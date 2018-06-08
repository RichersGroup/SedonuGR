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
	if(n_emit_core_per_bin>0)  emit_inner_source_by_bin();
	if(n_emit_zones_per_bin>0) emit_zones_by_bin();

	// sanity checks
	for(unsigned i=0; i<particles.size(); i++){
		if(particles.fate[i]==moving){
			for(unsigned j=0; j<4; j++){
				PRINT_ASSERT(particles.xup[j][i],==,particles.xup[j][i]);
				PRINT_ASSERT(particles.kup[j][i],==,particles.kup[j][i]);
			}
			PRINT_ASSERT(particles.N[i],==,particles.N[i]);
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

	const unsigned ns = species_list.size();
	const unsigned ng = grid->nu_grid_axis.size();
	const unsigned n_emit = ns*ng*n_emit_core_per_bin;
	unsigned n_emit_this_rank = n_emit / MPI_nprocs;
	if((int)(n_emit % MPI_nprocs) > MPI_myID) n_emit_this_rank++;
	particles.resize(size_before + n_emit_this_rank);

	unsigned n_created = 0;
	#pragma omp parallel for schedule(guided) collapse(3) reduction(+:n_created)
	for(unsigned s=0; s<ns; s++){
		for(unsigned g=0; g<ng; g++){
			for(int k=0; k<n_emit_core_per_bin; k++){
				unsigned global_id = k + n_emit_core_per_bin*g + n_emit_core_per_bin*ng*s;
				if((int)(global_id%MPI_nprocs) == MPI_myID){
					unsigned local_index = size_before + global_id/MPI_nprocs;
					create_surface_particle(particles, weight,s,g, local_index);
					if(particles.fate[local_index] == moving) n_created++;
				}
			}
		}
	}

	if(verbose) cout << "#   emit_inner_source_by_bin() created = " << n_created << " particles on rank 0 ("
			<< n_emit_this_rank-n_created << " rouletted during emission)" << endl;
}


//--------------------------------------------------------------------------
// emit particles due to viscous heating
//--------------------------------------------------------------------------
void Transport::emit_zones_by_bin(){
	int size_before = particles.size();
	double weight = 1./((double)n_emit_zones_per_bin);

	const unsigned ns = species_list.size();
	const unsigned ng = grid->nu_grid_axis.size();
	const unsigned nz = grid->rho.size();
	const unsigned n_emit = ns*ng*nz*n_emit_zones_per_bin;
	unsigned n_emit_this_rank = n_emit / MPI_nprocs;
	if((int)(n_emit % MPI_nprocs) > MPI_myID) n_emit_this_rank++;
	particles.resize(size_before + n_emit_this_rank);

	unsigned n_created = 0;
	#pragma omp parallel for reduction(+:n_created) schedule(guided) collapse(4)
	for (unsigned z_ind=0; z_ind<nz; z_ind++){
		for(unsigned s=0; s<ns; s++){
			for(unsigned g=0; g<ng; g++){
				for(int k=0; k<n_emit_zones_per_bin; k++){

					unsigned global_id = k + n_emit_zones_per_bin*g + n_emit_zones_per_bin*ng*s + n_emit_zones_per_bin*ng*ns*z_ind;
					if((int)(global_id%MPI_nprocs) == MPI_myID){
						unsigned local_index = size_before + global_id/MPI_nprocs;
						create_thermal_particle(particles, z_ind,weight,s,g, local_index);
						if(particles.fate[local_index] == moving){
							n_created++;
							for(unsigned d=0; d<4; d++) PRINT_ASSERT(particles.xup[d][local_index],==,particles.xup[d][local_index]);
						}
					}

				}
			}
		}
	}
	  
	double total_neutrinos = 0;
	for(unsigned i=0; i<species_list.size(); i++) total_neutrinos += N_net_emit[i];
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
void Transport::create_thermal_particle(ParticleList& output, const int z_ind,const double weight, const unsigned s, const unsigned g, const unsigned list_index)
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
		output.kup[3][list_index] = 0;
		output.N[list_index] = 0;
		output.fate[list_index] = rouletted;
		return;
	}

	// sample the frequency
	double nu3_min = pow(grid->nu_grid_axis.bottom(g), 3);
	double nu3_max = pow(grid->nu_grid_axis.top[g],    3);
	double nu3 = rangen.uniform(nu3_min, nu3_max);
	double nu = pow(nu3, 1./3.);

	// emit isotropically in comoving frame
	Tuple<double,4> kup_tet;
	isotropic_kup_tet(nu,kup_tet,&rangen);
	eh.set_kup_tet(kup_tet);
	update_eh_k_opac(&eh);

	// set the particle number
	eh.N = grid->BB[s][eh.eas_ind]/*.interpolate(eh.icube_spec)*/ * eh.absopac; // #/s/cm^3/sr/(Hz^3/3)
	eh.N *= grid->zone_lab_3volume(eh.z_ind);
	if(DO_GR) eh.N *= sqrt(eh.g.gammalow.det()) * (-eh.g.ndot(eh.u)); // comoving volume (d3x * volfac * Lorentz factor)
	eh.N *= weight * 1/*s*/ * 4.*pc::pi/*sr*/ * grid->nu_grid_axis.delta3(g)/3.0/*Hz^3/3*/;
	PRINT_ASSERT(eh.N,>=,0);
	PRINT_ASSERT(eh.N,<,1e99);
	eh.N0 = eh.N;

	// add to particle vector
	window(&eh);
	if(eh.fate == moving){
		PRINT_ASSERT(eh.N,>,0);

		// count up the emitted energy in each zone
		N_net_emit[eh.s] += eh.N;
		grid->l_emit[z_ind] -= eh.N * species_list[eh.s]->lepton_number;
		for(unsigned i=0; i<4; i++){
			grid->fourforce_emit[z_ind][i] -= eh.N * kup_tet[i];
		}
	}
	eh.get_Particle(output, list_index);
}


//------------------------------------------------------------
// General function to create a particle on the surface
// emitted isotropically outward in the comoving frame. 
//------------------------------------------------------------
bool reject_direction(const EinsteinHelper* eh, ThreadRNG* rangen){
	double xdotx = eh->g.dot<3>(eh->xup, eh->xup);
	double kdotk = eh->g.dot<3>(eh->kup, eh->kup);
	double xdotk = eh->g.dot<3>(eh->xup, eh->kup);
	double costheta = xdotk / sqrt(xdotx * kdotk);
	return (rangen->uniform() > costheta);
}
void Transport::create_surface_particle(ParticleList& output, const double weight, const unsigned int s, const unsigned int g, const unsigned list_index)
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
	double nu3_min = pow(grid->nu_grid_axis.bottom(g), 3);
	double nu3_max = pow(grid->nu_grid_axis.top[g],    3);
	double nu3 = rangen.uniform(nu3_min, nu3_max);
	double nu = pow(nu3, 1./3.);

	// sample outward direction
	Tuple<double,4> kup_tet;
	do{
		isotropic_kup_tet(nu,kup_tet,&rangen);
		eh.set_kup_tet(kup_tet);
	} while(reject_direction(&eh, &rangen));
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
			* weight;                    // 1/number of samples
	eh.N0 = eh.N;



	// add to particle vector
	window(&eh);
	if(eh.fate == moving){
		N_core_emit[eh.s] += eh.N;
	}
	eh.get_Particle(output, list_index);
}
