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
void Transport::which_event(const EinsteinHelper *eh, ParticleEvent *event, double* ds_com) const{
	PRINT_ASSERT(eh->N, >, 0);
	PRINT_ASSERT(eh->z_ind,>=,0);
	*event = nothing;

	// FIND D_ZONE= ====================================================================
	double d_zone = grid->zone_min_length(eh->z_ind) / sqrt(Metric::dot_Minkowski<3>(eh->kup,eh->kup)) * eh->kup_tet[3];
	d_zone = min(max(d_zone, d_zone*min_step_size), d_zone*max_step_size);
	PRINT_ASSERT(d_zone, >, 0);

	// FIND D_BOUNDARY
	double d_boundary = grid->d_boundary(*eh) * (1.0+TINY);
	d_boundary = max(d_boundary, d_zone*(1.0+TINY));
	*ds_com = min(d_boundary, d_zone);
	PRINT_ASSERT(d_boundary, >, 0);

	// FIND D_ABS
	if(eh->absopac*eh->ds_com > absorption_depth_limiter)
		*ds_com = absorption_depth_limiter/eh->absopac;

	// FIND D_RANDOMWALK
	double d_randomwalk = INFINITY;
	if(do_randomwalk && eh->scatopac*eh->ds_com>randomwalk_min_optical_depth){ // coarse check
		double D = pc::c / (3. * eh->scatopac);
		d_randomwalk = min(max(grid->d_randomwalk(*eh), d_zone*min_step_size), d_zone*max_step_size);
		if(r_core>0){
			// get a null test vector
			Tuple<double,4> ktest = -eh->xup;
			ktest[3] = 0;
			eh->g.normalize_null_changeupt(ktest);

			// limit d_randomwalk expecting movement towards core
			double r = radius(eh->xup);
			double kr = r;
			double kup_tet_t = -eh->g.dot<4>(ktest,eh->u);
			double ur = Metric::dot_Minkowski<3>(ktest,eh->u)/r;
			d_randomwalk = min(d_randomwalk, R_randomwalk(kr/kup_tet_t, ur, r-r_core, D));
		}
		if(eh->absopac > 0){
			double R_abs_limited = sqrt(absorption_depth_limiter * D / eh->absopac); // assumes a Delta t = 1s.
			d_randomwalk = min(d_randomwalk, R_abs_limited);
		}
		if(d_randomwalk == INFINITY) d_randomwalk = 1.1*randomwalk_min_optical_depth / eh->scatopac;
		PRINT_ASSERT(d_randomwalk,>=,0);
		if(eh->scatopac * d_randomwalk > randomwalk_min_optical_depth){ // real check
			*ds_com = d_randomwalk;
			*event = randomwalk;
		}
	}

	// FIND D_ELASTIC_SCATTER =================================================================
	double d_interact = INFINITY;
	if(*event!=randomwalk && eh->scatopac>0){
		double tau;
		do{
			tau = -log(rangen.uniform());
		} while(tau >= INFINITY);
		d_interact = tau / eh->scatopac;
		if(d_interact < *ds_com){
			*ds_com = d_interact;
			*event = elastic_scatter;
		}
	}

	// FIND D_INELASTIC_SCATTER =================================================================
	double d_inelastic_scatter = INFINITY;
	if(*event!=randomwalk && eh->inelastic_scatopac>0){
		double tau;
		do{
			tau = -log(rangen.uniform());
		} while(tau >= INFINITY);
		d_inelastic_scatter = tau / eh->inelastic_scatopac;
		if(d_inelastic_scatter < *ds_com){
			*ds_com = d_inelastic_scatter;
			*event = inelastic_scatter;
		}
	}
	PRINT_ASSERT(*ds_com, >=, 0);
	PRINT_ASSERT(*ds_com, <, INFINITY);
}

void Transport::move(EinsteinHelper *eh, bool do_absorption) const{
	PRINT_ASSERT(eh->ds_com,>=,0);
	PRINT_ASSERT(eh->N,>,0);
	PRINT_ASSERT(abs(eh->g.dot<4>(eh->kup,eh->kup)) / (eh->kup[3]*eh->kup[3]), <=, TINY);

	// save old values
	const EinsteinHelper eh_old = *eh;

	// convert ds_com into dlambda
	double dlambda = eh_old.ds_com / eh_old.kup_tet[3];
	PRINT_ASSERT(dlambda,>=,0);

	// get 2nd order x, 1st order estimate for k
	Tuple<double,4> order1 = eh_old.kup * dlambda;
	eh->xup += order1;
	if(DO_GR){
		Tuple<double,4> dk_dlambda = grid->dk_dlambda(eh_old);
		eh->kup = eh_old.kup + dk_dlambda * dlambda;
		Tuple<double,4> order2 = dk_dlambda * dlambda*dlambda * 0.5;
		if(radius(order2) < radius(order1)) eh->xup += order2;
	}

	// get new background data
	update_eh_background(eh);
	if(eh->fate==moving)
		update_eh_k_opac(eh);


	double tau=0, dN=0;
	if(do_absorption){

		// appropriately reduce the particle's energy from absorption
		// assumes kup_tet and absopac vary linearly along the trajectory
		tau = eh_old.ds_com * eh_old.absopac;
		eh->N *= exp(-tau);
		dN = eh_old.N - eh->N;
		window(eh);

		// store absorbed energy rate in *comoving* frame
		grid->fourforce_abs[eh_old.z_ind] += eh_old.kup_tet * dN/eh_old.zone_fourvolume;

		// store absorbed lepton number (same in both frames, except for the
		// factor of this_d which is divided out later
		if(species_list[eh->s]->lepton_number != 0){
			grid->l_abs[eh_old.z_ind] += dN * species_list[eh->s]->lepton_number / eh_old.zone_fourvolume;
		}
	}

	// tally in contribution to zone's distribution function (lab frame)
	// use old coordinates/directions to avoid problems with boundaries
	double avg_N = (tau>TINY ? dN/tau : (eh->N+eh_old.N)/2.);
	grid->distribution[eh_old.s]->count_single(eh_old.kup_tet, eh_old.dir_ind, avg_N*eh_old.ds_com*eh_old.kup_tet[3] / (eh_old.zone_fourvolume*pc::c));
	
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
		double ds_com;
		which_event(eh,&event, &ds_com);
		eh->ds_com = ds_com;
		PRINT_ASSERT(eh->ds_com ,>, 0);
		PRINT_ASSERT(eh->N,>,0);
		if(event==randomwalk)
		  random_walk(eh);
		else{
		  move(eh);
		  if(eh->z_ind>=0 and (event==elastic_scatter or event==inelastic_scatter))
		    scatter(eh, event);
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
		grid->spectrum[eh->s].count_single(eh->kup_tet, eh->dir_ind, e);
	}
	else if(eh->fate==absorbed)
		particle_core_abs_energy += e;
	else if(eh->fate==rouletted)
		particle_rouletted_energy += e;
	else assert(0);
}
