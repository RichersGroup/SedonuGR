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
	const unsigned nparticles = particles.size();

	//--- MOVE THE PARTICLES AROUND ---
	#pragma omp parallel for schedule(dynamic)
	for(unsigned i=0; i<nparticles; i++){
		// propagate each particle with an EinsteinHelper
		EinsteinHelper eh;
		eh.p = particles[i];
		update_eh_background(&eh);
		if(eh.z_ind >= 0) grid->interpolate_opacity(&eh);

		// save final data
		#pragma omp atomic
		n_active[eh.p.s]++;
		if(eh.p.fate == moving) propagate(&eh);
		if(eh.p.fate == escaped){
			const double nu = eh.nu(); // uses last-known metric
			double D[3] = {eh.p.kup[0], eh.p.kup[1], eh.p.kup[2]};
			Metric::normalize_Minkowski<3>(D);
			#pragma omp atomic
			n_escape[eh.p.s]++;
			#pragma omp atomic
			L_net_esc[eh.p.s] += eh.p.N * nu*pc::h;
			#pragma omp atomic
			N_net_esc[eh.p.s] += eh.p.N;
			grid->spectrum[eh.p.s].count(&eh, eh.p.N * nu*pc::h);
		}

		#pragma omp atomic
		ndone++;
		if(verbose && ndone%10000==0)
			cout << "\r"<<ndone<<"/"<<nparticles << " (" << (double)ndone/(double)nparticles*100<<"\%)";

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
	double d_zone = grid->zone_min_length(eh->z_ind); // REPLACE WITH detg^1/6
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

	if(r_core>0 && grid->radius(eh->p.xup)<r_core) eh->p.fate = absorbed;
	else if(eh->z_ind<0){
		grid->symmetry_boundaries(eh);
		update_eh_background(eh);
		if(eh->z_ind < 0) eh->p.fate = escaped;
	}
}

void Transport::tally_radiation(const EinsteinHelper *eh) const{
	PRINT_ASSERT(eh->z_ind, >=, 0);
	PRINT_ASSERT(eh->z_ind, <, grid->rho.size());
	PRINT_ASSERT(eh->ds_com, >=, 0);
	PRINT_ASSERT(eh->p.N, >, 0);
	PRINT_ASSERT(eh->nu(), >, 0);
	double to_add = 0;
	double decay_factor = 1.0 - exp(-eh->absopac * eh->ds_com); //same in both frames

	// tally in contribution to zone's distribution function (lab frame)
	if(eh->absopac>0) to_add = eh->p.N / eh->absopac * decay_factor;
	else to_add = eh->p.N * eh->ds_com;
	to_add *= eh->nu()*pc::h;
	PRINT_ASSERT(to_add,<,INFINITY);

	grid->distribution[eh->p.s]->count(eh, to_add);

	// store absorbed energy in *comoving* frame (will turn into rate by dividing by dt later)
	to_add = eh->p.N * decay_factor;
	PRINT_ASSERT(to_add,>=,0);

	Tuple<double,4> tmp_fourforce;
	for(unsigned i=0; i<4; i++){
		#pragma omp atomic
		grid->fourforce_abs[eh->z_ind][i] += eh->kup_tet[i] * to_add;
	}

	// store absorbed lepton number (same in both frames, except for the
	// factor of this_d which is divided out later
	if(species_list[eh->p.s]->lepton_number != 0){
		to_add *= species_list[eh->p.s]->lepton_number;
		#pragma omp atomic
		grid->l_abs[eh->z_ind] += to_add;
	}
}

void Transport::move(EinsteinHelper *eh) const{
	PRINT_ASSERT(eh->ds_com,>=,0);
	// FOR SCHWARZSCHILD PATH TEST
	//for(unsigned i=0; i<4; i++) cout << eh->p.xup[i] << "\t";
	//for(unsigned i=0; i<4; i++) cout << eh->p.kup[i] << "\t";
	//for(unsigned i=0; i<4; i++) cout << eh->kup_tet[i] << "\t";
	//cout << eh->nu() << endl;

	// translate the particle
	eh->integrate_geodesic();

	// appropriately reduce the particle's energy
	eh->p.N *= exp(-eh->absopac * eh->ds_com);
	window(eh);
	if(eh->p.fate==moving) PRINT_ASSERT(eh->p.N,>,0);

	update_eh_background(eh);
	if(eh->z_ind >= 0) grid->interpolate_opacity(eh);
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
		PRINT_ASSERT(eh->kup_tet[3],<,INFINITY);
		for(unsigned i=0; i<NDIMS; i++) PRINT_ASSERT(eh->dir_ind[i],<,grid->rho.axes[i].size());

		// decide which event happens
		which_event(eh,&event);

		// move particle the distance
		PRINT_ASSERT(eh->p.N,>,0);
		switch(event){
		case nothing:
			tally_radiation(eh);
			move(eh);
			break;
		case randomwalk:
			random_walk(eh);
			break;
		case interact:
			tally_radiation(eh);
			move(eh);
			if(eh->z_ind>0) scatter(eh);
			break;
		default:
			assert(0);
		}

		if(eh->p.fate==moving) window(eh);
		if(eh->p.fate==moving) boundary_conditions(eh);
		PRINT_ASSERT(eh->p.N,<,1e99);
	}

	// copy particle back out of LorentzHelper
	PRINT_ASSERT(eh->p.fate, !=, moving);
}
