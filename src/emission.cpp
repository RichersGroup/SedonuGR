#include "global_options.h"
#include <omp.h>
#include <math.h>
#include <gsl/gsl_sf_fermi_dirac.h>
#include "thread_RNG.h"
#include "transport.h"
#include "species_general.h"

//------------------------------------------------------------
// emit new particles
//------------------------------------------------------------
void transport::emit_particles(const double lab_dt)
{
	assert(lab_dt > 0);

	// complain if we're out of room for particles
	assert(n_emit_core  >= 0 || n_emit_zones >= 0);
	int n_emit = n_emit_core + n_emit_zones;
	if (total_particles() + n_emit > max_particles){
		if(rank0){
			cout << "Total particles: " << total_particles() << endl;
			cout << "n_emit: " << n_emit << endl;
			cout << "max_particles: " << max_particles << endl;
			cout << "# ERROR: Not enough particle space\n";
		}
		exit(10);
	}

	// set how many partices are going to bin-wise emission
	int n_emit_core_per_bin  = ratio_emit_by_bin * n_emit_core  / (number_of_bins()               );
	int n_emit_zones_per_bin = ratio_emit_by_bin * n_emit_zones / (number_of_bins()*grid->z.size());
	if(verbose && rank0 && ratio_emit_by_bin>0){
		cout << "#   " << n_emit_core_per_bin  << " particles per core bin" << endl;
		cout << "#   " << n_emit_zones_per_bin << " particles per zone bin" << endl;
	}


	// emit from the core and/or the zones
	if(verbose && rank0) cout << "# Emitting particles...";
	if(n_emit_core >0) emit_inner_source(n_emit_core_per_bin, lab_dt);
	if(n_emit_zones>0) emit_zones(n_emit_zones_per_bin, lab_dt);
	if(verbose && rank0) cout << "finished." << endl;
}

//------------------------------------------------------------
// begin with a blackbody distribution of particles
//------------------------------------------------------------
void transport::initialize_blackbody(double T, double munue){
	// remove existing particles
	particles.resize(0);

	for(unsigned z_ind=0; z_ind<grid->z.size(); z_ind++){
		// prep the emissivities
		for(unsigned s=0; s<species_list.size(); s++)
			species_list[s]->set_emis_to_BB_edens(T,munue*(double)species_list[s]->lepton_number);

		emit_zones(n_initial,1.0,t_now); // this is a hack - 1.0 causes perfect BB energy density
		grid->z[z_ind].e_emit = 0;
		grid->z[z_ind].l_emit = 0;
	}
}

//------------------------------------------------------------
// inject particles from a central luminous source
// Currently written to emit photons with 
// blackblody spectrum based on T_core and L_core
//------------------------------------------------------------
void transport::emit_inner_source(const int n_emit_core_per_bin, const double lab_dt, double t)
{
	assert(lab_dt > 0);
	assert(core_species_luminosity.N > 0);
	int n_emit_core_binwise = n_emit_core_per_bin*number_of_bins();
	int n_emit_core_regular = n_emit_core - n_emit_core_binwise;
	assert(n_emit_core_regular>=0);

	if(n_emit_core_regular>0 && ratio_emit_by_bin<1.0){
		const double Ep  = core_species_luminosity.N * lab_dt/n_emit_core_regular * (1.0 - ratio_emit_by_bin);
		assert(Ep > 0);
		assert(r_core>0);

        #pragma omp parallel for
		for (int i=0; i<n_emit_core_regular; i++){
			if(t<0) t = t_now + lab_dt*rangen.uniform();
			create_surface_particle(Ep,t);
		}
	}

	if(n_emit_core_per_bin > 0){
        #pragma omp parallel for
		for(unsigned s=0; s<species_list.size(); s++){
			for(unsigned g=0; g<species_list[s]->core_emis.size(); g++){
				double emis = species_list[s]->core_emis.get_value(g)*species_list[s]->core_emis.N;
				double bin_Ep = emis*lab_dt/n_emit_core_per_bin * ratio_emit_by_bin;
				if(bin_Ep>0) for(int k=0; k<n_emit_core_per_bin; k++){
					if(t<0) t = t_now + lab_dt*rangen.uniform();
					create_surface_particle(bin_Ep,t,s,g);
				}
			}
		}
	}
}


//--------------------------------------------------------------------------
// emit particles due to viscous heating
//--------------------------------------------------------------------------
void transport::emit_zones(const int n_emit_zones_per_bin, const double lab_dt, double t){
	assert(lab_dt > 0);

	unsigned gridsize = grid->z.size();
	double tmp_net_energy = 0;
	int n_emit_zones_regular = n_emit_zones - n_emit_zones_per_bin*grid->z.size()*number_of_bins();
	assert(n_emit_zones_regular>=0);


	// at this point therm means either viscous heating or regular emission, according to the logic above
    //#pragma omp parallel
	{
		// determine the net luminosity of each emission type over the whole grid
		// proper normalization due to frames not physical - just for distributing particles somewhat reasonably
	    #pragma omp parallel for reduction(+:tmp_net_energy)
		for(unsigned z_ind=0; z_ind<gridsize; z_ind++){
			double com_emit_energy = zone_comoving_therm_emit_energy(z_ind,lab_dt);
			tmp_net_energy += com_emit_energy;
		}
		assert(tmp_net_energy>0);

		// actually emit the particles in each zone
        //#pragma omp for schedule(guided)
		for (unsigned z_ind=0; z_ind<gridsize; z_ind++)
		{
			double com_emit_energy = zone_comoving_therm_emit_energy(z_ind,lab_dt);

			if(n_emit_zones_regular>0 && ratio_emit_by_bin<1.0){
				// how much this zone emits. Always emits correct energy even if number of particles doesn't add up.
				unsigned this_n_emit = (double)n_emit_zones_regular * (com_emit_energy / tmp_net_energy) + 0.5;
				if(com_emit_energy>0 && this_n_emit==0) this_n_emit = 1;
				double Ep = com_emit_energy / (double)this_n_emit * (1.0 - ratio_emit_by_bin);
				for (unsigned k=0; k<this_n_emit; k++){
					if(t<0) t = t_now + lab_dt*rangen.uniform();
					create_thermal_particle(z_ind,Ep,t);
				}
			}


			// create particles in all bins
			if(n_emit_zones_per_bin > 0){
				for(unsigned s=0; s<species_list.size(); s++){
					for(unsigned g=0; g<species_list[s]->number_of_bins(); g++){
						double bin_Ep = bin_comoving_therm_emit_energy(z_ind, s, g, lab_dt)/n_emit_zones_per_bin * ratio_emit_by_bin;
						if(bin_Ep>0) for(int k=0; k<n_emit_zones_per_bin; k++){
							if(t<0) t = t_now + lab_dt*rangen.uniform();
							create_thermal_particle(z_ind,bin_Ep,t,s,g);
						}
					}
				}
			}

			// record emissivity
			grid->z[z_ind].e_emit += com_emit_energy;
			grid->z[z_ind].l_emit += zone_comoving_therm_emit_leptons(z_ind,lab_dt);
		}// loop over zones
	}// #pragma omp parallel
}



//----------------------------------------------------------------------------------------
// Helper functions for emit_zones
//----------------------------------------------------------------------------------------

// return the cell's luminosity from thermal emission (erg/s, comoving frame)
double transport::zone_comoving_therm_emit_energy(const int z_ind, const double lab_dt) const{
	if(!grid->good_zone(z_ind)) return 0;                 // don't emit from superluminal zones
	else if(grid->zone_radius(z_ind) < r_core) return 0; // don't emit within core
	else{
		double H=0;
		double four_vol = grid->zone_lab_volume(z_ind) * lab_dt; //relativistic invariant - same in comoving frame.
		for(unsigned i=0; i<species_list.size(); i++){
			double species_lum = species_list[i]->integrate_zone_emis(z_ind) * 4*pc::pi * four_vol;
			assert(species_lum >= 0);
			H += species_lum;
		}
		assert(H>=0);
		return H;
	}
}

double transport::bin_comoving_therm_emit_energy(const int z_ind, const int s, const int g, const double lab_dt) const{
	if(!grid->good_zone(z_ind)) return 0;                 // don't emit from superluminal zones
	else if(grid->zone_radius(z_ind) < r_core) return 0; // don't emit within core
	else{
		double four_vol = grid->zone_lab_volume(z_ind) * lab_dt; //relativistic invariant - same in comoving frame.
		double H = species_list[s]->bin_emis(z_ind,g) * 4*pc::pi * four_vol;
		assert(H>=0);
		return H;
	}

}

// return the cell's luminosity from thermal emission (erg/s, comoving frame)
double transport::zone_comoving_therm_emit_leptons(const int z_ind, const double lab_dt) const{
	if(!grid->good_zone(z_ind)) return 0; //don't emit from superluminal zones
	else{
		double L=0;
		double four_vol = grid->zone_lab_volume(z_ind) * lab_dt; //relativistic invariant - same in comoving frame.
		for(unsigned i=0; i<species_list.size(); i++){
			double species_lum = species_list[i]->integrate_zone_lepton_emis(z_ind) * 4*pc::pi * four_vol;
			L += species_lum;
		}
		return L;
	}
}

//------------------------------------------------------------
// General function to create a particle in zone i
// emitted isotropically in the comoving frame. 
// Useful for thermal radiation emitted all througout
// the grid
//------------------------------------------------------------
void transport::create_thermal_particle(const int z_ind, const double Ep, const double t, const int s, const int g)
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
	normalize(p.D);

	// sample the species and frequency
	p.s = (s>=0 ? s : sample_zone_species(z_ind) );
	assert(p.s >= 0);
	assert(p.s < (int)species_list.size());
	p.nu = species_list[p.s]->sample_zone_nu(z_ind,g);
	assert(p.nu > 0);

	// lorentz transform from the comoving to lab frame
	transform_comoving_to_lab(&p,z_ind);

	// add to particle vector
	assert(particles.size() < particles.capacity());
    //#pragma omp critical
	particles.push_back(p);

	// count up the emitted energy in each zone
    #pragma omp atomic
	L_net_lab[p.s] += p.e;
	#pragma omp atomic
	E_avg_lab[p.s] += p.nu * p.e;
}


//------------------------------------------------------------
// General function to create a particle on the surface
// emitted isotropically outward in the comoving frame. 
//------------------------------------------------------------
void transport::create_surface_particle(const double Ep, const double t, const int s, const int g)
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
	normalize(p.D);

	// get index of current zone
	assert(grid->zone_index(p.x) >= 0);

	// sample the species and frequency
	p.s = (s>=0 ? s : sample_core_species());
	assert(p.s >= 0);
	assert(p.s < (int)species_list.size());
	p.nu = species_list[p.s]->sample_core_nu(g);

	// lorentz transform from the comoving to lab frame
	transform_comoving_to_lab(&p);

	// add to particle vector
    #pragma omp critical
	particles.push_back(p);
    //#pragma omp atomic
	//L_net_lab[p.s] += p.e;
}
