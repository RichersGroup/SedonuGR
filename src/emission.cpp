/*
//  Copyright (c) 2015, California Institute of Technology and the Regents
//  of the University of California, based on research sponsored by the
//  United States Department of Energy. All rights reserved.
//
//  This file is part of Sedonu.
//
//  Sedonu is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  Neither the name of the California Institute of Technology (Caltech)
//  nor the University of California nor the names of its contributors 
//  may be used to endorse or promote products derived from this software
//  without specific prior written permission.
//
//  Sedonu is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with Sedonu.  If not, see <http://www.gnu.org/licenses/>.
//
*/

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
void transport::emit_particles()
{
	// complain if we're out of room for particles
	assert(n_emit_core  >= 0 || n_emit_zones >= 0 || n_emit_core_per_bin>0 || n_emit_zones_per_bin>=0);
	int n_emit = do_emit_by_bin ?
			(n_emit_core_per_bin + n_emit_zones_per_bin*grid->z.size()) * species_list.size()*number_of_bins() :
			n_emit_core + n_emit_zones;
	if (total_particles() + n_emit > max_particles){
		if(rank0){
			cout << "Total particles: " << total_particles() << endl;
			cout << "n_emit: " << n_emit << endl;
			cout << "max_particles: " << max_particles << endl;
			cout << "# ERROR: Not enough particle space\n";
		}
		exit(10);
	}


	// emit from the core and/or the zones
	if(verbose && rank0) cout << "# Emitting particles..." << endl;
	if(do_emit_by_bin){
		if(n_emit_core_per_bin>0)  emit_inner_source_by_bin();
		if(n_emit_zones_per_bin>0) emit_zones_by_bin();
	}
	else{
		if(n_emit_core>0)  emit_inner_source();
		if(n_emit_zones>0) emit_zones();
	}
}

//------------------------------------------------------------
// inject particles from a central luminous source
// Currently written to emit photons with 
// blackblody spectrum based on T_core and L_core
//------------------------------------------------------------
void transport::emit_inner_source_by_bin(){
	int size_before = particles.size();

	double avgEp = 0;
	#pragma omp parallel for reduction(+:avgEp)
	for(unsigned s=0; s<species_list.size(); s++){
		for(unsigned g=0; g<species_list[s]->core_emis.size(); g++){
			double emis = species_list[s]->core_emis.get_value(g)*species_list[s]->core_emis.N;
			double bin_Ep = emis/n_emit_core_per_bin * do_emit_by_bin; // assume lab_dt = 1.0
			if(bin_Ep>0){
				for(int k=0; k<n_emit_core_per_bin; k++) create_surface_particle(bin_Ep,s,g);
				avgEp += bin_Ep * n_emit_core_per_bin;
			}
		}
	}
	int n_created = particles.size()-size_before;
	avgEp /= (double)n_created;

	PRINT_ASSERT(n_created,==,particles.size()-size_before);
	if(verbose && rank0) cout << "#   <E_p> (emit_inner_source_by_bin) = " << avgEp << " erg ("
			<< n_created << " particles)" << endl;
}
void transport::emit_inner_source()
{
	PRINT_ASSERT(core_species_luminosity.N,>,0);
	int size_before = particles.size();

	const double Ep  = core_species_luminosity.N / (double)n_emit_core; // assume lab_dt = 1.0
	PRINT_ASSERT(Ep,>,0);
	PRINT_ASSERT(r_core,>,0);

	#pragma omp parallel for
	for (int i=0; i<n_emit_core; i++) create_surface_particle(Ep);

	if(verbose && rank0) cout << "#   <E_p> (emit_inner_source) = " << Ep << " erg ("
			<< particles.size()-size_before << " particles)" << endl;
}


//--------------------------------------------------------------------------
// emit particles due to viscous heating
//--------------------------------------------------------------------------
void transport::emit_zones_by_bin(){
	double avgEp = 0;
	int size_before = particles.size();

	#pragma omp parallel for reduction(+:avgEp)
	for (unsigned z_ind=0; z_ind<grid->z.size(); z_ind++) if(grid->zone_radius(z_ind) >= r_core){

		double com_emit_energy = zone_comoving_therm_emit_energy(z_ind);

		for(unsigned s=0; s<species_list.size(); s++){
			for(unsigned g=0; g<species_list[s]->number_of_bins(); g++){
				double bin_Ep = bin_comoving_therm_emit_energy(z_ind, s, g)/n_emit_zones_per_bin;
				if(bin_Ep>0){
					for(int k=0; k<n_emit_zones_per_bin; k++) create_thermal_particle(z_ind,bin_Ep,s,g);
					avgEp += bin_Ep * n_emit_zones_per_bin;
				}
			}
		}

		// record emissivity
		#pragma omp atomic
		grid->z[z_ind].e_emit += com_emit_energy;
		#pragma omp atomic
		grid->z[z_ind].l_emit += zone_comoving_therm_emit_leptons(z_ind);
	}

	int n_created = particles.size() - size_before;
	avgEp /= (double)n_created;
	if(verbose && rank0) cout << "#   <E_p> (emit_zones_by_bin) = " << avgEp << " erg ("
			<< particles.size()-size_before << " particles)" << endl;
}
void transport::emit_zones(){
	double tmp_net_energy = 0;
	int size_before = particles.size();
	double avgEp = 0;


	// determine the net luminosity of each emission type over the whole grid
	// proper normalization due to frames not physical - just for distributing particles somewhat reasonably
	#pragma omp parallel for reduction(+:tmp_net_energy)
	for(unsigned z_ind=0; z_ind<grid->z.size(); z_ind++){
		double com_biased_emit_energy = zone_comoving_biased_therm_emit_energy(z_ind);
		tmp_net_energy += com_biased_emit_energy;
	}
	PRINT_ASSERT(tmp_net_energy,>,0);

	// actually emit the particles in each zone
	#pragma omp parallel for schedule(guided) reduction(+:avgEp)
	for (unsigned z_ind=0; z_ind<grid->z.size(); z_ind++) if(grid->zone_radius(z_ind) >= r_core)
	{
		double com_biased_emit_energy = zone_comoving_biased_therm_emit_energy(z_ind);
		double com_emit_energy = zone_comoving_biased_therm_emit_energy(z_ind);

		// how much this zone emits. Always emits correct energy even if number of particles doesn't add up.
		// randomly emit additional particle based on remainder of allocation after emitting an integer number of packets
		if(com_biased_emit_energy>0){
			double this_zones_share = (double)n_emit_zones * (com_biased_emit_energy / tmp_net_energy);
			unsigned this_n_emit = (unsigned)this_zones_share;
			double remainder = this_zones_share - (double)this_n_emit;
			PRINT_ASSERT(remainder,<,1.0);
			PRINT_ASSERT(remainder,>=,0.0);
			if(rangen.uniform()<remainder) this_n_emit += 1;

			double Ep = com_emit_energy / (double)this_n_emit;
			avgEp += com_emit_energy;
			for (unsigned k=0; k<this_n_emit; k++) create_thermal_particle(z_ind,Ep);
		}

		// record emissivity
		#pragma omp atomic
		grid->z[z_ind].e_emit += com_biased_emit_energy;
		#pragma omp atomic
		grid->z[z_ind].l_emit += zone_comoving_therm_emit_leptons(z_ind);
	}// loop over zones

	int n_created = particles.size() - size_before;
	avgEp /= (double)n_created;
	if(verbose && rank0) cout << "#   <E_p> (emit_zones) = " << avgEp << " erg ("
			<< n_created << " particles)" << endl;
}



//----------------------------------------------------------------------------------------
// Helper functions for emit_zones
//----------------------------------------------------------------------------------------

// return the cell's luminosity from thermal emission (erg/s, comoving frame)
double transport::zone_comoving_therm_emit_energy(const int z_ind) const{
	if(!grid->good_zone(z_ind)) return 0;                 // don't emit from superluminal zones
	else if(grid->zone_radius(z_ind) < r_core) return 0; // don't emit within core
	else{
		double H=0;
		double four_vol = grid->zone_lab_volume(z_ind); //relativistic invariant - same in comoving frame. Assume lab_dt=1.0
		for(unsigned i=0; i<species_list.size(); i++){
			double species_lum = species_list[i]->integrate_zone_emis(z_ind) * 4*pc::pi * four_vol;
			PRINT_ASSERT(species_lum,>=,0);
			H += species_lum;
		}
		PRINT_ASSERT(H,>=,0);
		return H;
	}
}
double transport::zone_comoving_biased_therm_emit_energy(const int z_ind) const{
	if(!grid->good_zone(z_ind)) return 0;                 // don't emit from superluminal zones
	else if(grid->zone_radius(z_ind) < r_core) return 0; // don't emit within core
	else{
		double H=0;
		double four_vol = grid->zone_lab_volume(z_ind); //relativistic invariant - same in comoving frame. Assume lab_dt=1.0
		for(unsigned i=0; i<species_list.size(); i++){
			double species_lum = species_list[i]->integrate_zone_biased_emis(z_ind) * 4*pc::pi * four_vol;
			PRINT_ASSERT(species_lum,>=,0);
			H += species_lum;
		}
		PRINT_ASSERT(H,>=,0);
		return H;
	}
}
double transport::bin_comoving_therm_emit_energy(const int z_ind, const int s, const int g) const{
	if(!grid->good_zone(z_ind)) return 0;                 // don't emit from superluminal zones
	else if(grid->zone_radius(z_ind) < r_core) return 0; // don't emit within core
	else{
		double four_vol = grid->zone_lab_volume(z_ind); //relativistic invariant - same in comoving frame. Assume lab_dt=1.0
		double H = species_list[s]->bin_emis(z_ind,g) * 4*pc::pi * four_vol;
		PRINT_ASSERT(H,>=,0);
		return H;
	}

}

// return the cell's luminosity from thermal emission (erg/s, comoving frame)
double transport::zone_comoving_therm_emit_leptons(const int z_ind) const{
	if(!grid->good_zone(z_ind)) return 0; //don't emit from superluminal zones
	else{
		double L=0;
		double four_vol = grid->zone_lab_volume(z_ind); //relativistic invariant - same in comoving frame. Assume lab_dt=1.0
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
void transport::create_thermal_particle(const int z_ind, const double Ep, const int s, const int g)
{
	PRINT_ASSERT(Ep,>,0);
	PRINT_ASSERT(z_ind,>=,0);
	PRINT_ASSERT(z_ind,<,(int)grid->z.size());

	// basic particle properties
	particle p;
	p.fate = moving;
	p.e  = Ep;

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
	if(s>=0) p.s = s;
	else sample_zone_species(&p,z_ind);
	PRINT_ASSERT(p.s,>=,0);
	PRINT_ASSERT(p.s,<,(int)species_list.size());
	species_list[p.s]->sample_zone_nu(p,z_ind,g);
	PRINT_ASSERT(p.nu,>,0);

	// lorentz transform from the comoving to lab frame
	transform_comoving_to_lab(&p,z_ind);

	// sample tau
	double lab_opac=0, abs_frac=0, dshift_l2c=0;
	lab_opacity(&p,z_ind,&lab_opac,&abs_frac,&dshift_l2c);
	sample_tau(&p,lab_opac,abs_frac);
	window(&p,z_ind);

	// add to particle vector
	if(p.fate == moving){
		PRINT_ASSERT(particles.size(),<,particles.capacity());
		PRINT_ASSERT(p.e,>,0);
		PRINT_ASSERT(p.tau,>,0);
		#pragma omp critical
		particles.push_back(p);

		// count up the emitted energy in each zone
		#pragma omp atomic
		L_net_lab[p.s] += p.e;
		#pragma omp atomic
		E_avg_lab[p.s] += p.nu * p.e;
		#pragma omp atomic
		N_net_lab[p.s] += p.e / (p.nu * pc::h);
	}
}


//------------------------------------------------------------
// General function to create a particle on the surface
// emitted isotropically outward in the comoving frame. 
//------------------------------------------------------------
void transport::create_surface_particle(const double Ep, const int s, const int g)
{
	PRINT_ASSERT(Ep,>,0);

	// set basic properties
	particle p;
	p.fate = moving;
	p.e = Ep;

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

	// get index of current zone
	const int z_ind = grid->zone_index(p.x);
	PRINT_ASSERT(z_ind,>=,0);

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


	// sample the species and frequency
	p.s = (s>=0 ? s : sample_core_species());
	PRINT_ASSERT(p.s,>=,0);
	PRINT_ASSERT(p.s,<,(int)species_list.size());
	p.nu = species_list[p.s]->sample_core_nu(g);

	// sample tau
	double lab_opac=0, abs_frac=0, dshift_l2c=0;
	lab_opacity(&p,z_ind,&lab_opac,&abs_frac,&dshift_l2c);
	sample_tau(&p,lab_opac,abs_frac);
	window(&p,z_ind);

	// add to particle vector
	if(p.fate == moving){
		PRINT_ASSERT(particles.size(),<,particles.capacity());
	    #pragma omp critical
		particles.push_back(p);
	    #pragma omp atomic
		L_core_lab[p.s] += p.e;
	}
}
