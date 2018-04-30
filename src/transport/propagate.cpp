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

	unsigned ndone=0;
	unsigned last_percent = 0;
	const unsigned nparticles = particles.size();

	//--- MOVE THE PARTICLES AROUND ---
	#pragma omp parallel for schedule(dynamic)
	for(unsigned i=0; i<nparticles; i++){
		// propagate each particle with an EinsteinHelper
		EinsteinHelper eh;
		eh.p = particles[i];
		eh.N0 = eh.p.N;
		update_eh_background(&eh);
		update_eh_k_opac(&eh);

		// save final data
		#pragma omp atomic
		n_active[eh.p.s]++;
		if(eh.p.fate == moving) propagate(&eh);
		if(eh.p.fate == escaped){
			const double nu = eh.nu(); // uses last-known metric
			#pragma omp atomic
			n_escape[eh.p.s]++;
			#pragma omp atomic
			L_net_esc[eh.p.s] += eh.p.N * nu*pc::h;
			#pragma omp atomic
			N_net_esc[eh.p.s] += eh.p.N;
			grid->spectrum[eh.p.s].count(eh.kup_tet, eh.dir_ind, eh.p.N * nu*pc::h);
		}

		if(verbose){
			#pragma omp atomic
			ndone++;
			unsigned this_percent = (double)ndone/(double)nparticles*100.;
			if(this_percent > last_percent){
				last_percent = this_percent;
				#pragma omp critical
				cout << "\r"<<ndone<<"/"<<nparticles << " (" << last_percent<<"\%)" << flush;
			}
		}
		PRINT_ASSERT(eh.p.fate, !=, moving);
		particles[i] = eh.p;
	} //#pragma omp parallel for
	cout << endl;


	double tot=0,core=0,roulette=0,esc=0;
	#pragma omp parallel for reduction(+:tot,core,roulette,esc)
	for(unsigned i=0; i<particles.size(); i++){
		PRINT_ASSERT(particles[i].fate,!=,moving);
		const double e  = particles[i].N * particles[i].kup[3];
		if(particles[i].fate!=rouletted) tot       += e;
		if(particles[i].fate==escaped  ) esc       += e;
		if(particles[i].fate==absorbed ) core      += e;
		if(particles[i].fate==rouletted) roulette  += e;
	}
	particle_total_energy = tot;
	particle_core_abs_energy = core;
	particle_rouletted_energy = roulette;
	particle_escape_energy = esc;

	// remove the dead particles
	particles.resize(0);
}

//--------------------------------------------------------
// Decide what happens to the particle
//--------------------------------------------------------
void Transport::which_event(EinsteinHelper *eh, ParticleEvent *event) const{
	PRINT_ASSERT(eh->p.N, >, 0);
	PRINT_ASSERT(eh->z_ind,>=,0);

	// FIND D_ZONE= ====================================================================
	double d_boundary = grid->d_boundary(eh);
	double d_zone = grid->zone_min_length(eh->z_ind) / sqrt(Metric::dot_Minkowski<3>(eh->p.kup,eh->p.kup)) * eh->kup_tet[3];
	d_boundary = min(max(d_boundary, d_zone*min_step_size), d_zone*max_step_size);
	PRINT_ASSERT(d_zone, >, 0);
	*event = nothing;
	eh->ds_com = d_boundary;

	// FIND D_RANDOMWALK
	double d_randomwalk = INFINITY;
	if(do_randomwalk && eh->scatopac*eh->ds_com>randomwalk_min_optical_depth){ // coarse check
		d_randomwalk = grid->d_randomwalk(eh);
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


void Transport::boundary_conditions(EinsteinHelper *eh) const{
	PRINT_ASSERT(eh->p.fate,==,moving);

	if(r_core>0 && radius(eh->p.xup)<r_core) eh->p.fate = absorbed;
	else if(eh->z_ind<0){
		grid->symmetry_boundaries(eh);
		update_eh_background(eh);
		update_eh_k_opac(eh);
		if(eh->z_ind < 0) eh->p.fate = escaped;
	}
}

void Transport::move(EinsteinHelper *eh) const{
	PRINT_ASSERT(eh->ds_com,>=,0);
	PRINT_ASSERT(abs(eh->g.dot<4>(eh->p.kup,eh->p.kup)) / (eh->p.kup[3]*eh->p.kup[3]), <=, TINY);
	// FOR SCHWARZSCHILD PATH TEST
	//for(unsigned i=0; i<4; i++) cout << eh->p.xup[i] << "\t";
	//for(unsigned i=0; i<4; i++) cout << eh->p.kup[i] << "\t";
	//for(unsigned i=0; i<4; i++) cout << eh->kup_tet[i] << "\t";
	//cout << eh->g.gtt << "\t";
	//cout << eh->nu() << endl;

	// save old values
	double old_kup[4], old_kup_tet[4];
	for(unsigned i=0; i<4; i++){
		old_kup[i] = eh->p.kup[i];
		old_kup_tet[i] = eh->kup_tet[i];
	}
	double old_absopac = eh->absopac;
	double old_N = eh->p.N;
	unsigned old_z_ind = eh->z_ind;
	unsigned old_dir_ind[NDIMS+1];
	for(unsigned i=0; i<NDIMS+1; i++) old_dir_ind[i] = eh->dir_ind[i];

	// convert ds_com into dlambda
	double dlambda = eh->ds_com / eh->kup_tet[3];
	PRINT_ASSERT(dlambda,>=,0);

	// get dk_dlambda at current position
	double dk_dlambda[4] = {0,0,0,0};
	if(DO_GR){
		eh->christoffel.contract2(eh->p.kup,dk_dlambda);
		for(unsigned i=0; i<4; i++) dk_dlambda[i] *= -1;
	}

	// get 2nd order x, 1st order estimate for k
	double kup1[4];
	for(unsigned i=0; i<4; i++){
		kup1[i] = old_kup[i] + dk_dlambda[i]*dlambda;
		double order1 = old_kup[i]*dlambda;
		double order2 = 0.5*dk_dlambda[i]*dlambda*dlambda;
		eh->p.xup[i] += order1 + (abs(order2/order1)<1. ? order2 : 0);
		PRINT_ASSERT(eh->p.xup[i],==,eh->p.xup[i]);
	}

	// get new background data
	update_eh_background(eh);

	// apply second order correction to k
	if(DO_GR and eh->z_ind>0){
		double dk_dlambda_2[4], kup2[4];
		eh->g.normalize_null_preserveupt(eh->p.kup);
		eh->christoffel.contract2(eh->p.kup,dk_dlambda_2);
		for(unsigned i=0; i<4; i++){
			dk_dlambda_2[i] *= -1;
			kup2[i] = old_kup[i] + dlambda*dk_dlambda_2[i];
			eh->p.kup[i] = 0.5 * (kup1[i] + kup2[i]);
			PRINT_ASSERT(eh->p.kup[i],==,eh->p.kup[i]);
		}
	}
	PRINT_ASSERT(abs(eh->p.kup[3]),<,INFINITY);
	update_eh_k_opac(eh);

	// appropriately reduce the particle's energy from absorption
	// assumes kup_tet and absopac vary linearly along the trajectory
	double ds_com_new = dlambda*eh->kup_tet[3];
	double tau1 = 1./3. * (eh->ds_com*old_absopac + ds_com_new*eh->absopac);
	double tau2 = 1./6. * (eh->ds_com*eh->absopac + ds_com_new*old_absopac);
	double tau = tau1 + tau2;
	eh->p.N *= exp(-tau);
	double dN = old_N - eh->p.N;
	window(eh);

	// get average k
	double avg_kup_tet[4];
	for(unsigned i=0; i<4; i++) avg_kup_tet[i] = 0.5 * (old_kup_tet[i] + eh->kup_tet[i]);

	// tally in contribution to zone's distribution function (lab frame)
	// use average N and integrate assuming kt_tet varies linearly over the step
	double avg_N = (tau>TINY ? dN/tau : old_N);
	double integral_k2dlambda = dlambda/3. * (old_kup_tet[3]*old_kup_tet[3] + old_kup_tet[3]*eh->kup_tet[3] + eh->kup_tet[3]*eh->kup_tet[3]);
	grid->distribution[eh->p.s]->count(avg_kup_tet, old_dir_ind, avg_N*integral_k2dlambda);

	// store absorbed energy in *comoving* frame (will turn into rate by dividing by dt later)
	for(unsigned i=0; i<4; i++){
		#pragma omp atomic
		grid->fourforce_abs[old_z_ind][i] += dN * avg_kup_tet[i];
	}

	// store absorbed lepton number (same in both frames, except for the
	// factor of this_d which is divided out later
	if(species_list[eh->p.s]->lepton_number != 0){
		#pragma omp atomic
		grid->l_abs[old_z_ind] += dN * species_list[eh->p.s]->lepton_number;
	}
}


//--------------------------------------------------------
// Propagate a single monte carlo particle until
// it  escapes, is absorbed, or the time step ends
//--------------------------------------------------------
void Transport::propagate(EinsteinHelper *eh) const{
	ParticleEvent event;

	PRINT_ASSERT(eh->p.fate, ==, moving);

	while (eh->p.fate == moving)
	{
		PRINT_ASSERT(eh->z_ind,>=,0);
		PRINT_ASSERT(eh->p.N,>,0);
		PRINT_ASSERT(eh->p.N,<,1e99);
		PRINT_ASSERT(eh->p.kup[3],>,0);
		PRINT_ASSERT(eh->kup_tet[3],>,0);
		PRINT_ASSERT(eh->p.kup[3],<,INFINITY);
		PRINT_ASSERT(eh->nu(),<,10.*grid->nu_grid_axis.max());
		for(unsigned i=0; i<NDIMS; i++) PRINT_ASSERT(eh->dir_ind[i],<,grid->rho.axes[i].size());

		// decide which event happens
		which_event(eh,&event);

		// move particle the distance
		PRINT_ASSERT(eh->p.N,>,0);
		switch(event){
		case nothing:
			move(eh);
			break;
		case randomwalk:
			random_walk(eh);
			break;
		case interact:
			move(eh);
			if(eh->z_ind>0) scatter(eh);
			break;
		default:
			assert(0);
		}

		if(eh->p.fate==moving) window(eh);
		if(eh->p.fate==moving) boundary_conditions(eh);
		if(eh->p.fate==moving) PRINT_ASSERT(abs(eh->g.dot<4>(eh->p.kup,eh->p.kup)) / (eh->p.kup[3]*eh->p.kup[3]), <=, TINY);

		PRINT_ASSERT(eh->p.N,<,1e99);
	}
}
