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
#include "Transport.h"
#include "Species.h"
#include "Grid.h"
#include <cstring>
#include "EinsteinHelper.h"

using namespace std;
namespace pc = physical_constants;

void Transport::propagate_particles()
{
	if(verbose) cout << "# Propagating particles..." << endl;

	size_t ndone=0;
	size_t last_percent = 0;
	const size_t nparticles = particles.size();

	//--- MOVE THE PARTICLES AROUND ---
	#pragma omp parallel for schedule(dynamic)
	for(size_t i=0; i<nparticles; i++){
		if(particles[i].fate != moving){
			#pragma omp atomic
			ndone++;
			continue;
		}

		// propagate each particle with an EinsteinHelper
		EinsteinHelper eh;
		eh.set_Particle(particles[i]);
		eh.N0 = eh.N;
		update_eh_background(&eh);
		update_eh_k_opac(&eh);
		propagate(&eh);

		if(verbose){
			#pragma omp atomic
			ndone++;
			size_t this_percent = (double)ndone/(double)nparticles*100.;
			if(this_percent > last_percent){
				last_percent = this_percent;
				#pragma omp critical
				cout << "\r"<<ndone<<"/"<<nparticles << " (" << last_percent<<"%)" << flush;
			}
		}
		PRINT_ASSERT(eh.fate, !=, moving);
		particles[i] = eh.get_Particle();
	} //#pragma omp parallel for
	if(verbose) cout << endl;

	// remove the dead particles, erase the memory
	particles.resize(0);
}

//--------------------------------------------------------
// Decide what happens to the particle
//--------------------------------------------------------
void Transport::which_event(EinsteinHelper *eh, ParticleEvent *event) const{
	PRINT_ASSERT(eh->N, >, 0);
	PRINT_ASSERT(eh->z_ind,>=,0);

	// FIND D_ZONE= ====================================================================
	double d_boundary = grid->d_boundary(*eh);
	double d_zone = grid->zone_min_length(eh->z_ind) / sqrt(Metric::dot_Minkowski<3>(eh->kup,eh->kup)) * eh->kup_tet[3];
	d_boundary = min(max(d_boundary, d_zone*min_step_size), d_zone*max_step_size);
	PRINT_ASSERT(d_zone, >, 0);
	*event = nothing;
	eh->ds_com = d_boundary;

	// FIND D_RANDOMWALK
	double d_randomwalk = INFINITY;
	if(do_randomwalk && eh->scatopac*eh->ds_com>randomwalk_min_optical_depth){ // coarse check
		d_randomwalk = grid->d_randomwalk(*eh);
		if(eh->absopac > 0){
			double D = pc::c / (3. * eh->scatopac);
			double R_abs_limited = sqrt(randomwalk_absorption_depth_limit * D / eh->absopac);
			d_randomwalk = min(d_randomwalk, R_abs_limited);
		}
		if(d_randomwalk == INFINITY) d_randomwalk = 1.1*randomwalk_min_optical_depth / eh->scatopac;
		PRINT_ASSERT(d_randomwalk,>=,0);
		if(eh->scatopac * d_randomwalk > randomwalk_min_optical_depth){ // real check
			eh->ds_com = d_randomwalk;
			*event = randomwalk;
		}
	}

	// FIND D_INTERACT =================================================================
	double d_interact = INFINITY;
	if(*event==nothing && eh->scatopac>0){
		double tau;
		do{
			tau = -log(rangen.uniform());
		} while(tau >= INFINITY);
		d_interact = tau / eh->scatopac;
		if(d_interact < eh->ds_com){
			eh->ds_com = d_interact;
			*event = interact;
		}
	}
	PRINT_ASSERT(eh->ds_com, >=, 0);
	PRINT_ASSERT(eh->ds_com, <, INFINITY);
}

void Transport::move(EinsteinHelper *eh, bool do_absorption) const{
	PRINT_ASSERT(eh->ds_com,>=,0);
	PRINT_ASSERT(eh->N,>,0);
	PRINT_ASSERT(abs(eh->g.dot<4>(eh->kup,eh->kup)) / (eh->kup[3]*eh->kup[3]), <=, TINY);

	// save old values
	Tuple<double,4> old_kup = eh->kup;
	Tuple<double,4> old_kup_tet = eh->kup_tet;
	double old_absopac = eh->absopac;
	double old_N = eh->N;
	size_t old_z_ind = eh->z_ind;
	size_t old_dir_ind[NDIMS+1];
	for(size_t i=0; i<NDIMS+1; i++) old_dir_ind[i] = eh->dir_ind[i];

	// convert ds_com into dlambda
	double dlambda = eh->ds_com / eh->kup_tet[3];
	PRINT_ASSERT(dlambda,>=,0);

	// get 2nd order x, 1st order estimate for k
	Tuple<double,4> order1 = old_kup * dlambda;
	for(size_t i=0; i<4; i++)
		eh->xup[i] += order1[i];
	if(DO_GR){
		Tuple<double,4> dk_dlambda = grid->dk_dlambda(*eh);
		Tuple<double,4> order2 = dk_dlambda * dlambda*dlambda * 0.5;
		eh->kup = old_kup + dk_dlambda * dlambda;
		for(size_t i=0; i<4; i++)
			eh->xup[i] += (abs(order2[i]/order1[i])<1. ? order2[i] : 0);
	}

	// get new background data
	update_eh_background(eh);
	if(eh->fate==moving)
		update_eh_k_opac(eh);


	double tau=0, dN=0;
	if(do_absorption){

		// appropriately reduce the particle's energy from absorption
		// assumes kup_tet and absopac vary linearly along the trajectory
		double ds_com_new = dlambda*eh->kup_tet[3];
		double tau1 = 1./3. * (eh->ds_com*old_absopac + ds_com_new*eh->absopac);
		double tau2 = 1./6. * (eh->ds_com*eh->absopac + ds_com_new*old_absopac);
		tau = tau1 + tau2;
		eh->N *= exp(-tau);
		dN = old_N - eh->N;
		window(eh);

		// store absorbed energy in *comoving* frame (will turn into rate by dividing by dt later)
		for(size_t i=0; i<4; i++){
			grid->fourforce_abs[old_z_ind][i] += dN * old_kup_tet[i];
		}

		// store absorbed lepton number (same in both frames, except for the
		// factor of this_d which is divided out later
		if(species_list[eh->s]->lepton_number != 0){
			grid->l_abs[old_z_ind] += dN * species_list[eh->s]->lepton_number;
		}
	}

	// tally in contribution to zone's distribution function (lab frame)
	// use old coordinates/directions to avoid problems with boundaries
	double avg_N = (tau>TINY ? dN/tau : old_N);
	grid->distribution[eh->s]->count(old_kup_tet, old_dir_ind, avg_N*dlambda*old_kup_tet[3]*old_kup_tet[3]);


}


//--------------------------------------------------------
// Propagate a single monte carlo particle until
// it  escapes, is absorbed, or the time step ends
//--------------------------------------------------------
void Transport::propagate(EinsteinHelper *eh){
	ParticleEvent event;

	PRINT_ASSERT(eh->fate, ==, moving);
	n_active[eh->s]++;

	while (eh->fate == moving)
	{
		PRINT_ASSERT(eh->z_ind,>=,0);
		PRINT_ASSERT(eh->N,>,0);
		PRINT_ASSERT(eh->N,<,1e99);
		PRINT_ASSERT(eh->kup[3],>,0);
		PRINT_ASSERT(eh->kup_tet[3],>,0);
		PRINT_ASSERT(eh->kup[3],<,INFINITY);
		for(size_t i=0; i<NDIMS; i++) PRINT_ASSERT(eh->dir_ind[i],<,grid->rho.axes[i].size());

		// decide which event happens
		which_event(eh,&event);

		// move particle the distance
		PRINT_ASSERT(eh->N,>,0);
		switch(event){
		case nothing:
			move(eh);
			break;
		case randomwalk:
			random_walk(eh);
			break;
		case interact:
			move(eh);
			if(eh->z_ind>0 && eh->scatopac>0) scatter(eh);
			break;
		default:
			assert(0);
		}

		if(eh->fate==moving) window(eh);
		if(eh->fate==moving) PRINT_ASSERT(abs(eh->g.dot<4>(eh->kup,eh->kup)) / (eh->kup[3]*eh->kup[3]), <=, TINY);

		PRINT_ASSERT(eh->N,<,1e99);
	}

	PRINT_ASSERT(eh->fate,!=,moving);
	double e = eh->N * eh->kup[3];
	if(eh->fate==escaped){
		PRINT_ASSERT(e,>=,0);
		particle_escape_energy += e;
		n_escape[eh->s]++;
		L_net_esc[eh->s] += e;
		N_net_esc[eh->s] += eh->N;
		grid->spectrum[eh->s].count(eh->kup_tet, eh->dir_ind, e);
	}
	else if(eh->fate==absorbed)
		particle_core_abs_energy += e;
	else if(eh->fate==rouletted)
		particle_rouletted_energy += e;
	else assert(0);
}
