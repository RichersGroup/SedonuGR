#include <omp.h>
#include <math.h>
#include <gsl/gsl_sf_fermi_dirac.h>
#include "thread_RNG.h"
#include "transport.h"
#include "species_general.h"
#include "global_options.h"

//------------------------------------------------------------
// emit new particles
//------------------------------------------------------------
void transport::emit_particles(const double dt)
{
	assert(dt > 0);

	// complain if we're out of room for particles
	assert(n_emit_core >= 0);
	assert(n_emit_therm >= 0);
	assert(n_emit_decay >= 0);
	int n_emit = n_emit_core + n_emit_therm + n_emit_decay;
	if (total_particles() + n_emit > max_particles){
		cout << "# ERROR: Not enough particle space\n";
		exit(10);
	}

	// emit from the core and/or the zones
	if(n_emit_core >0) emit_inner_source(dt);
	if(n_emit_therm>0) emit_zones(dt, n_emit_therm, &transport::zone_therm_lum,      &transport::create_thermal_particle);
	if(n_emit_decay>0) emit_zones(dt, n_emit_decay, &transport::zone_decay_lum,      &transport::create_decay_particle);
}


//------------------------------------------------------------
// inject particles from a central luminous source
// Currently written to emit photons with 
// blackblody spectrum based on T_core and L_core
//------------------------------------------------------------
void transport::emit_inner_source(const double dt)
{
	assert(dt > 0);
	assert(core_species_luminosity.N > 0);
	const double Ep  = core_species_luminosity.N * dt/n_emit_core;
	L_net += core_species_luminosity.N;
	assert(Ep > 0);

#pragma omp parallel for
	for (int i=0; i<n_emit_core; i++){
		double t = t_now + dt*rangen.uniform();
		create_surface_particle(Ep,t);
	}
}


//--------------------------------------------------------------------------
// emit particles due to viscous heating
//--------------------------------------------------------------------------
void transport::emit_zones(const double dt,
		const int n_emit,
		double (transport::*zone_lum)(const int) const,
		void (transport::*create_particle)(const int,const double,const double)){

	assert(dt > 0);
	assert(n_emit >= 0);

	unsigned gridsize = grid->z.size();
	double tmp_net_lum = NaN;
	double Ep = NaN;

	// at this point therm means either viscous heating or regular emission, according to the logic above
#pragma omp parallel
	{
		// determine the net luminosity of each emission type over the whole grid
#pragma omp for reduction(+:tmp_net_lum)
		for(unsigned i=0; i<gridsize; i++){
			double zone_luminosity = (this->*zone_lum)(i);
			assert(zone_luminosity >= 0);
			tmp_net_lum += zone_luminosity;
		}

		// determine the energy of each emitted particle
#pragma omp single
		{
			Ep = tmp_net_lum*dt / (double)n_emit;
			L_net += tmp_net_lum;
		}

		// actually emit the particles in each zone
#pragma omp for schedule(guided)
		for (unsigned i=0; i<gridsize; i++)
		{
			// this zone's luminosity and number of emitted particles.
			// randomly decide whether last particle gets added based on the remainder.
			double almost_n_emit  = (double)n_emit * (this->*zone_lum)(i)/tmp_net_lum;
			unsigned this_n_emit = (int)almost_n_emit + (int)( rangen.uniform() < fmod(almost_n_emit,1.0) );
			assert(this_n_emit >= 0);

			// create the particles
			double t;
			for (unsigned k=0; k<this_n_emit; k++){
				t = t_now + dt*rangen.uniform();
				(this->*create_particle)(i,Ep,t);
			}
		}// loop over zones
	}// #pragma omp parallel
}



//----------------------------------------------------------------------------------------
// Helper functions for emit_zones
//----------------------------------------------------------------------------------------

// return the cell's luminosity from thermal emission (erg/s)
double transport::zone_therm_lum(const int zone_index) const{
	double H=0;
	for(unsigned i=0; i<species_list.size(); i++){
		double species_lum = species_list[i]->integrate_zone_emis(zone_index) * 4*pc::pi * grid->zone_volume(zone_index);
		assert(species_lum >= 0);
		H += species_lum;
	}
	return H;
}

// return the total decay luminosity from the cell (erg/s)
double transport::zone_decay_lum(const int zone_index) const{
	if(rank0) cout << "WARNING: emission in zones from decay is not yet implemented!" << endl;
	return 0;
}


//------------------------------------------------------------
// General function to create a particle in zone i
// emitted isotropically in the comoving frame. 
// Useful for thermal radiation emitted all througout
// the grid
//------------------------------------------------------------
void transport::create_thermal_particle(const int z_ind, const double Ep, const double t)
{
	assert(Ep > 0);
	assert(z_ind >= 0);
	assert(z_ind < (int)grid->z.size());

	// basic particle properties
	particle p;
	p.e  = Ep;
	p.t  = t;

	// random sample position in zone
	vector<double> rand(3,0);
	rand[0] = rangen.uniform();
	rand[1] = rangen.uniform();
	rand[2] = rangen.uniform();
	vector<double> r;
	grid->cartesian_sample_in_zone(z_ind,rand,r);
	p.x[0] = r[0];
	p.x[1] = r[1];
	p.x[2] = r[2];

	// emit isotropically in comoving frame
	double mu  = 1 - 2.0*rangen.uniform();
	double phi = 2.0*pc::pi*rangen.uniform();
	double smu = sqrt(1 - mu*mu);
	p.D[0] = smu*cos(phi);
	p.D[1] = smu*sin(phi);
	p.D[2] = mu;
	p.normalize_direction();

	// sample the species and frequency
	int s = sample_zone_species(z_ind);
	p.s = s;
	p.nu = species_list[s]->sample_zone_nu(z_ind);
	assert(p.nu > 0);

	// lorentz transform from the comoving to lab frame
	transform_comoving_to_lab(&p);

	// add to particle vector
#pragma omp critical
	particles.push_back(p);

	// count up the emitted energy/leptons in each zone
#pragma omp atomic
	grid->z[z_ind].e_emit += p.e;
#pragma omp atomic
	grid->z[z_ind].l_emit += p.e/(pc::h*p.nu) * (double)species_list[s]->lepton_number;
}


//------------------------------------------------------------
// General function to create a particle on the surface
// emitted isotropically outward in the comoving frame. 
//------------------------------------------------------------
void transport::create_surface_particle(const double Ep, const double t)
{
	assert(Ep > 0);

	// set basic properties
	particle p;
	p.e = Ep;
	p.t = t;

	// pick initial position on photosphere
	double phi_core   = 2*pc::pi*rangen.uniform();
	double cosp_core  = cos(phi_core);
	double sinp_core  = sin(phi_core);
	double cost_core  = 1 - 2.0*rangen.uniform();
	double sint_core  = sqrt(1-cost_core*cost_core);
	// double spatial coordinates
	double a_phot = r_core + r_core*1e-10;
	p.x[0] = a_phot*sint_core*cosp_core;
	p.x[1] = a_phot*sint_core*sinp_core;
	p.x[2] = a_phot*cost_core;

	// pick photon propagation direction wtr to local normal
	double phi_loc = 2*pc::pi*rangen.uniform();
	// choose sqrt(R) to get outward, cos(theta) emission
	double cost_loc  = sqrt(rangen.uniform());
	double sint_loc  = sqrt(1 - cost_loc*cost_loc);
	// local direction vector
	double D_xl = sint_loc*cos(phi_loc);
	double D_yl = sint_loc*sin(phi_loc);
	double D_zl = cost_loc;
	// apply rotation matrix to convert D vector into overall frame
	p.D[0] = cost_core*cosp_core*D_xl-sinp_core*D_yl+sint_core*cosp_core*D_zl;
	p.D[1] = cost_core*sinp_core*D_xl+cosp_core*D_yl+sint_core*sinp_core*D_zl;
	p.D[2] = -sint_core*D_xl+cost_core*D_zl;
	p.normalize_direction();

	// get index of current zone
	assert(grid->zone_index(p.x) >= 0);

	// sample the species and frequency
	int s = sample_core_species();
	assert(s >= 0);
	assert(s < (int)species_list.size());
	p.s = s;
	p.nu = species_list[s]->sample_core_nu();

	// lorentz transform from the comoving to lab frame
	transform_comoving_to_lab(&p);

	// add to particle vector
#pragma omp critical
	particles.push_back(p);
}


//------------------------------------------------------------
// General function to create a particle from radioactive decay
//------------------------------------------------------------
void transport::create_decay_particle(const int zone_index, const double Ep, const double t)
{
	if(rank0) cout << "WARNING: transport::create_decay_particle is not yet implemented!" << endl;
}
