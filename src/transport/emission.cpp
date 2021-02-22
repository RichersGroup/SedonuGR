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
	// emit from the core and/or the zones
	if(verbose) cout << "# Emitting particles..." << endl;
	if(n_emit_core_per_bin>0 and r_core>0)  emit_inner_source_by_bin();
	if(n_emit_zones_per_bin>0) emit_zones_by_bin();

	// sanity checks
	for(size_t i=0; i<particles.size(); i++){
		if(particles[i].fate==moving){
			for(size_t j=0; j<4; j++){
				PRINT_ASSERT(particles[i].xup[j],==,particles[i].xup[j]);
				PRINT_ASSERT(particles[i].kup[j],==,particles[i].kup[j]);
			}
			PRINT_ASSERT(particles[i].N,==,particles[i].N);
		}
	}
}

//------------------------------------------------------------
// inject particles from a central luminous source
// Currently written to emit photons with 
// blackblody spectrum based on T_core and L_core
//------------------------------------------------------------
void Transport::emit_inner_source_by_bin(){
	int size_before = particles.size();
	double weight = 1./((double)n_emit_core_per_bin);

	const size_t ns = species_list.size();
	const size_t ng = grid->nu_grid_axis.size();
	const size_t n_emit = ns*ng*n_emit_core_per_bin;
	size_t n_emit_this_rank = n_emit / MPI_nprocs;
	if((int)(n_emit % MPI_nprocs) > MPI_myID) n_emit_this_rank++;
	particles.resize(size_before + n_emit_this_rank);

	size_t n_created = 0;
	#pragma omp parallel for schedule(guided) collapse(3) reduction(+:n_created)
	for(size_t s=0; s<ns; s++){
		for(size_t g=0; g<ng; g++){
			for(int k=0; k<n_emit_core_per_bin; k++){
				size_t global_id = k + n_emit_core_per_bin*g + n_emit_core_per_bin*ng*s;
				if((int)(global_id%MPI_nprocs) == MPI_myID){
					size_t local_index = size_before + global_id/MPI_nprocs;
					particles[local_index] = create_surface_particle(weight,s,g);
					if(particles[local_index].fate == moving) n_created++;
				}
			}
		}
	}

	if(verbose) cout << "#   emit_inner_source_by_bin() created = " << n_created << " particles on rank 0 ("
			<< n_emit_this_rank-n_created << " rouletted during emission)" << endl;
}


//--------------------------------------------------------------------------
// emit particles
//--------------------------------------------------------------------------
void Transport::emit_zones_by_bin(){
	int size_before = particles.size();
	double weight = 1./((double)n_emit_zones_per_bin);

	const size_t ns = species_list.size();
	const size_t ng = grid->nu_grid_axis.size();
	const size_t nz = grid->rho.size();
	const size_t n_emit = ns*ng*nz*n_emit_zones_per_bin;
	size_t n_emit_this_rank = n_emit / MPI_nprocs;
	if((int)(n_emit % MPI_nprocs) > MPI_myID) n_emit_this_rank++;
	particles.resize(size_before + n_emit_this_rank);

	size_t n_created = 0;
	#pragma omp parallel for reduction(+:n_created) schedule(guided) collapse(4)
	for (size_t z_ind=0; z_ind<nz; z_ind++){
		for(size_t s=0; s<ns; s++){
			for(size_t g=0; g<ng; g++){
				for(int k=0; k<n_emit_zones_per_bin; k++){

					size_t global_id = k + n_emit_zones_per_bin*g + n_emit_zones_per_bin*ng*s + n_emit_zones_per_bin*ng*ns*z_ind;
					if((int)(global_id%MPI_nprocs) == MPI_myID){
						size_t local_index = size_before + global_id/MPI_nprocs;
						particles[local_index] = create_thermal_particle(z_ind,weight,s,g);
						if(particles[local_index].fate == moving){
							n_created++;
							for(size_t d=0; d<4; d++) PRINT_ASSERT(particles[local_index].xup[d],==,particles[local_index].xup[d]);
						}
					}

				}
			}
		}
	}
	  
	double total_neutrinos = 0;
	for(size_t i=0; i<species_list.size(); i++) total_neutrinos += N_net_emit[i];
	if(verbose) cout << "#   emit_zones_by_bin() created " << n_created << " particles on rank 0 ("
			<< total_neutrinos << " neutrinos) ("
			<< n_emit_this_rank-n_created << " rouletted immediately)" << endl;
}


//------------------------------------------------------------
// General function to create a particle in zone i
// emitted isotropically in the comoving frame. 
// Useful for thermal radiation emitted all througout
// the grid
//------------------------------------------------------------
Particle Transport::create_thermal_particle(const int z_ind,const double weight, const size_t s, const size_t g)
{
	PRINT_ASSERT(z_ind,>=,0);
	PRINT_ASSERT(z_ind,<,(int)grid->rho.size());
	PRINT_ASSERT(s,<,species_list.size());
	
	EinsteinHelper eh;
	eh.fate = moving;
	eh.s = s;

	// random sample position in zone
	eh.xup = grid->sample_in_zone(z_ind,&rangen);
	eh.xup[3] = 0;
	update_eh_background(&eh);
	if(eh.z_ind<0 || radius(eh.xup)<r_core){
		Particle output;
		output.kup[3] = 0;
		output.N = 0;
		output.fate = rouletted;
		return output;
	}

	// sample the frequency
	double nu=0;
	while(nu==0){ // reject nu=0
		double nu3 = rangen.uniform( pow(grid->nu_grid_axis.bottom(g),3), pow(grid->nu_grid_axis.top[g],3) );
		nu = pow(nu3, 1./3.);
	}
	PRINT_ASSERT(nu,>,0);
	nu = grid->nu_grid_axis.mid[g];

	// emit isotropically in comoving frame
	Tuple<double,4> kup_tet;
	kup_tet[3] = nu * pc::h;
	isotropic_kup_tet(kup_tet,&rangen);
	eh.set_kup_tet(kup_tet);
	update_eh_k_opac(&eh);

	// set the particle number
	double T = grid->T.interpolate(eh.icube_vol);
	double mu = grid->munue.interpolate(eh.icube_vol) * species_list[s]->lepton_number;
	eh.N = number_blackbody(T,mu,nu) * eh.absopac * species_list[s]->weight; // #/s/cm^3/sr/(Hz^3/3)
	eh.N *= eh.zone_fourvolume;// frame-independent four-volume
	eh.N *= weight * 4.*pc::pi/*sr*/ * grid->nu_grid_axis.delta3(g)/3.0/*Hz^3/3*/;
	PRINT_ASSERT(eh.N,>=,0);
	PRINT_ASSERT(eh.N,<,1e99);
	eh.N0 = eh.N;

	// add to particle vector
	window(&eh);
	if(eh.fate == moving){
		PRINT_ASSERT(eh.N,>,0);

		// count up the emitted energy in each zone
		N_net_emit[eh.s] += eh.N;
		grid->l_emit[z_ind] -= eh.N * species_list[eh.s]->lepton_number / eh.zone_fourvolume;
		for(size_t i=0; i<4; i++){
			grid->fourforce_emit[z_ind][i] -= eh.N * kup_tet[i] / eh.zone_fourvolume;
		}
	}
	return eh.get_Particle();
}


//------------------------------------------------------------
// General function to create a particle on the surface
// emitted isotropically outward in the comoving frame. 
//------------------------------------------------------------
Particle Transport::create_surface_particle(const double weight, const size_t s, const size_t g)
{
	PRINT_ASSERT(weight,>,0);
	PRINT_ASSERT(weight,!=,INFINITY);
	PRINT_ASSERT(s,<,species_list.size());

	EinsteinHelper eh;
	eh.fate = moving;
	eh.s = s;

	// pick initial position on photosphere
	random_core_x(eh.xup);
	eh.xup[3] = 0;
	update_eh_background(&eh);

	// sample the frequency
	double nu=0;
	while(nu==0){ // reject nu=0
		double nu3 = rangen.uniform( pow(grid->nu_grid_axis.bottom(g),3), pow(grid->nu_grid_axis.top[g],3) );
		nu = pow(nu3, 1./3.);
	}
	PRINT_ASSERT(nu,>,0);

	// sample outward direction
	Tuple<double,4> kup_tet;
	kup_tet[3] = nu * pc::h;
	double costheta;
	do{
		isotropic_kup_tet(kup_tet,&rangen);
		eh.set_kup_tet(kup_tet);
		costheta = eh.g.dot<3>(eh.xup, eh.kup) / sqrt(eh.g.dot<3>(eh.kup, eh.kup) * eh.g.dot<3>(eh.xup, eh.xup));
	} while(r_core>0 and reject_direction(costheta, 2.)); // 2. makes pdf = costheta
	update_eh_k_opac(&eh);

	//get the number of neutrinos in the particle
	double T = species_list[s]->T_core;
	double mu = species_list[s]->mu_core;
	double multiplier = species_list[s]->core_lum_multiplier * species_list[s]->weight;
	PRINT_ASSERT(nu,>,0);
	eh.N = number_blackbody(T,mu,nu)   // #/s/cm^2/sr/(Hz^3/3)
			* 1                          //   s
			* (4.0*pc::pi*r_core*r_core) //     cm^2
			* pc::pi                     //          sr (including factor of 1/2 for integrating over cos(theta)
			* grid->nu_grid_axis.delta3(g)/3.0 //        Hz^3/3
			* multiplier                 // overall scaling
			* (DO_GR ? eh.g.alpha : 1.)  // time lapse (s)
			* weight;                    // 1/number of samples
	eh.N0 = eh.N;

	// add to particle vector
	window(&eh);
	if(eh.fate == moving){
		N_core_emit[eh.s] += eh.N;
	}
	return eh.get_Particle();
}
