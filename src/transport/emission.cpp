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

#include <omp.h>
#include "Transport.h"
#include "Species.h"
#include "Grid.h"
#include "global_options.h"

using namespace std;
namespace pc = physical_constants;

//------------------------------------------------------------
// emit new particles
//------------------------------------------------------------
void Transport::emit_particles()
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
void Transport::emit_inner_source_by_bin(){
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

	PRINT_ASSERT(n_created,==,(int)particles.size()-size_before);
	if(verbose && rank0) cout << "#   <E_p> (emit_inner_source_by_bin) = " << avgEp << " erg ("
			<< n_created << " particles)" << endl;
}
void Transport::emit_inner_source()
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
void Transport::emit_zones_by_bin(){
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
void Transport::emit_zones(){
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
			unsigned this_n_emit = (unsigned)this_zones_share; // truncate.
			double remainder = this_zones_share - (double)this_n_emit;
			PRINT_ASSERT(remainder,<,1.0);
			PRINT_ASSERT(remainder,>=,0.0);

			if(rangen.uniform()<remainder) this_n_emit += 1;
			double Ep = com_emit_energy / (double)this_n_emit;
			if(this_n_emit>0){
				PRINT_ASSERT(Ep,>,0);
				avgEp += Ep * (double)this_n_emit;

				// keep emitted energy statistically constant
				if(this_zones_share<1){
					PRINT_ASSERT(this_n_emit,==,1);
					Ep /= this_zones_share;
				}
			}

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
double Transport::zone_comoving_therm_emit_energy(const int z_ind) const{
	if(z_ind<0) return 0;
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
double Transport::zone_comoving_biased_therm_emit_energy(const int z_ind) const{
	if(z_ind<0) return 0;
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
double Transport::bin_comoving_therm_emit_energy(const int z_ind, const int s, const int g) const{
	if(z_ind<0) return 0;
	else if(grid->zone_radius(z_ind) < r_core) return 0; // don't emit within core
	else{
		double four_vol = grid->zone_lab_volume(z_ind); //relativistic invariant - same in comoving frame. Assume lab_dt=1.0
		double H = species_list[s]->bin_emis(z_ind,g) * 4*pc::pi * four_vol;
		PRINT_ASSERT(H,>=,0);
		return H;
	}

}

// return the cell's luminosity from thermal emission (erg/s, comoving frame)
double Transport::zone_comoving_therm_emit_leptons(const int z_ind) const{
	if(z_ind<0) return 0;
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
void Transport::create_thermal_particle(const int z_ind, const double Ep, const int s, const int g)
{
	PRINT_ASSERT(Ep,>,0);
	PRINT_ASSERT(z_ind,>=,0);
	PRINT_ASSERT(z_ind,<,(int)grid->z.size());

	Particle pcom;
	pcom.fate = moving;
	pcom.e = Ep;

	// random sample position in zone
	double rand[3];
	rand[0] = rangen.uniform();
	rand[1] = rangen.uniform();
	rand[2] = rangen.uniform();
	grid->sample_in_zone(z_ind,rand,3,pcom.xup,3);
	pcom.xup[3] = 0;

	// get the velocity
	double v[3];
	grid->interpolate_fluid_velocity(pcom.xup,3,v,3,z_ind);

	// species
	pcom.s = s>=0 ? s : sample_zone_species(z_ind,&pcom.e);
	PRINT_ASSERT(pcom.s,>=,0);
	PRINT_ASSERT(pcom.s,<,(int)species_list.size());

	// frequency
	double nu = species_list[pcom.s]->sample_zone_nu(z_ind,&pcom.e,g);

	// emit isotropically in comoving frame
	grid->isotropic_kup(nu,pcom.kup,pcom.xup,4,&rangen);

	// set up LorentzHelper
	LorentzHelper lh(exponential_decay);
	lh.set_v(v,3);
	lh.set_p<com>(&pcom);

	// sample tau
	get_opacity(&lh,z_ind);
	sample_tau(&lh);
	window(&lh,z_ind);

	// add to particle vector
	if(lh.p_fate() == moving){
		PRINT_ASSERT(particles.size(),<,particles.capacity());
		PRINT_ASSERT(lh.p_e(lab),>,0);
		PRINT_ASSERT(lh.p_tau(),>,0);
		#pragma omp critical
		particles.push_back(lh.particle_copy(lab));

		// count up the emitted energy in each zone
		#pragma omp atomic
		L_net_lab[lh.p_s()] += lh.p_e(lab);
		#pragma omp atomic
		N_net_lab[lh.p_s()] += lh.p_e(lab) / (lh.p_nu(lab) * pc::h);
	}
}


//------------------------------------------------------------
// General function to create a particle on the surface
// emitted isotropically outward in the comoving frame. 
//------------------------------------------------------------
void Transport::create_surface_particle(const double Ep, const int s, const int g)
{
	PRINT_ASSERT(Ep,>,0);

	Particle plab;
	plab.fate = moving;
	plab.e = Ep;

	// pick initial position on photosphere
	double D[3];
	grid->random_core_x_D(r_core,&rangen,plab.xup,D,3);
	plab.xup[3] = 0;

	// get index of current zone
	const int z_ind = grid->zone_index(plab.xup, 3);
	PRINT_ASSERT(z_ind,>=,0);

	// sample the species
	plab.s = ( s>=0 ? s : sample_core_species() );
	PRINT_ASSERT(plab.s,>=,0);
	PRINT_ASSERT(plab.s,<,(int)species_list.size());

	// sample the frequency
	plab.kup[3] = species_list[plab.s]->sample_core_nu(g) / pc::c * 2.0*pc::pi;
	for(int i=0; i<3; i++) plab.kup[i] = D[i] * plab.kup[3];

	// set up LorentzHelper
	LorentzHelper lh(exponential_decay);
	double v[3];
	grid->interpolate_fluid_velocity(plab.xup,3,v,3,z_ind);
	lh.set_v(v,3);
	lh.set_p<lab>(&plab);

	// sample the optical depth
	get_opacity(&lh,z_ind);
	sample_tau(&lh);
	window(&lh,z_ind);

	// add to particle vector
	if(lh.p_fate() == moving){
		PRINT_ASSERT(particles.size(),<,particles.capacity());
	    #pragma omp critical
		particles.push_back(lh.particle_copy(lab));
	    #pragma omp atomic
		L_core_lab[lh.p_s()] += lh.p_e(lab);
	}
}
