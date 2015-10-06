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
#include <stdio.h>
#include <math.h>
#include <iostream>
#include "transport.h"
#include "species_general.h"

void transport::propagate_particles(const double lab_dt)
{
	if(verbose && rank0) cout << "# Propagating particles...";

	//--- MOVE THE PARTICLES AROUND ---
    #pragma omp parallel for schedule(guided)
	for(unsigned i=0; i<particles.size(); i++){
		particle* p = &particles[i];
        #pragma omp atomic
		n_active[p->s]++;
		propagate(p,lab_dt);
		if(p->fate == escaped){
			#pragma omp atomic
			n_escape[p->s]++;
			#pragma omp atomic
			L_net_esc[p->s] += p->e;
			#pragma omp atomic
			E_avg_esc[p->s] += p->nu * p->e;
			#pragma omp atomic
			N_net_esc[p->s] += p->e / (p->nu*pc::h);
			species_list[p->s]->spectrum.count(p->D, p->nu, p->e);
		}
	} //#pragma omp parallel for

	double tot=0,core=0,fluid=0,esc=0;
	#pragma omp parallel for reduction(+:tot,core,fluid,esc)
	for(unsigned i=0; i<particles.size(); i++){
		tot += particles[i].e;
		if(particles[i].fate==escaped) esc += particles[i].e;
		if(particles[i].fate==absorbed){
			if(grid->zone_index(particles[i].x)>=0 ) fluid += particles[i].e;
			else core += particles[i].e;
		}
	}
	particle_total_energy += tot;
	particle_core_abs_energy += core;
	particle_fluid_abs_energy += fluid;
	particle_escape_energy += esc;

	// remove the dead particles
	remove_dead_particles();

	if(verbose && rank0) cout << "finished." << endl;
}

void transport::remove_dead_particles(){
	if(iterative) particles.resize(0);
	else{
		vector<particle>::iterator pIter = particles.begin();
		while(pIter != particles.end()){
			if(pIter->fate==absorbed || pIter->fate==escaped){
				*pIter = particles[particles.size()-1];
				particles.pop_back();
			}
			else pIter++;
		}
	}
}


//--------------------------------------------------------
// Decide what happens to the particle
//--------------------------------------------------------
void transport::which_event(const particle *p, const int z_ind, const double dt, const double lab_opac,
		double *d_smallest, ParticleEvent *event) const{
	assert(lab_opac >= 0);

	double d_zone     = numeric_limits<double>::infinity();
	double d_time     = numeric_limits<double>::infinity();
	double d_interact = numeric_limits<double>::infinity();
	double d_boundary = numeric_limits<double>::infinity();
	double tau_r = NaN;                               // random optical depth
	assert(z_ind >= -1);

	if(z_ind >= 0){ //i.e. within the simulation region
		// FIND D_ZONE= ====================================================================
		// maximum step size inside zone
		d_zone = step_size * grid->zone_min_length(z_ind);
		assert(d_zone > 0);

		// FIND D_INTERACT =================================================================
		// random optical depth to next interaction
		tau_r = -1.0*log(1.0 - rangen.uniform());
		// step size to next interaction event
		d_interact  = tau_r/lab_opac;
		if (lab_opac == 0) d_interact = numeric_limits<double>::infinity();
		assert(d_interact>=0);
	}

	// FIND D_TIME ====================================================================
	// find distance to end of time step
	if(dt>0){
		double tstop = t_now + dt;
		d_time = (tstop - p->t)*pc::c;
		assert(d_time > 0);
	}

	// FIND D_BOUNDARY ================================================================
	d_boundary = grid->lab_dist_to_boundary(p);
	assert(d_boundary > 0);

	// find out what event happens (shortest distance) =====================================
	*d_smallest = numeric_limits<double>::infinity();
	*event = interact;
	*d_smallest = d_interact;
	if( d_zone <= *d_smallest ){
		*event  = zoneEdge;
		*d_smallest = d_zone;
	}
	if( d_time <= *d_smallest ){
		*event  = timeStep;
		*d_smallest = d_time;
	}
	if( d_boundary <= *d_smallest ){
		*event = boundary;
		*d_smallest = d_boundary;
		if(z_ind >= 0) *d_smallest *= (1.0 + grid_general::tiny); // bump just over the boundary if in simulation domain
		else           *d_smallest *= (1.0 - grid_general::tiny); // don't overshoot outward through the inner boundary
	}
	assert(*d_smallest >= 0);
}


void transport::event_boundary(particle* p, const int z_ind) const{
	assert(z_ind==-1 || z_ind==-2);

	// if outside the domain
	if(z_ind == -2){
		int new_ind = z_ind;
		if(reflect_outer){
			grid->reflect_outer(p);
			assert(p->fate == moving);
			new_ind = grid->zone_index(p->x);
			assert(new_ind >= 0);
			assert(new_ind < (int)grid->z.size());
			assert(p->nu > 0);
		}
		else{
			grid->reflect_symmetry(p);
			new_ind = grid->zone_index(p->x);
		}

		if(new_ind < 0) p->fate = escaped;
	}

	// if inside the inner boundary
	if(z_ind==-1){
		if(p->r() < r_core) p->fate = absorbed;
		else if(p->x_dot_d() >= 0){
			// set the particle just outside the inner boundary
			cout << "ERROR: have not yet implemented passing out through the inner boundary without overshooting" << endl;
			exit(5);
		}
		else assert(p->fate == moving); // the particle just went into the inner boundary
	}
}

void transport::tally_radiation(const particle* p, const int z_ind, const double dshift_l2c, const double lab_d, const double lab_opac, const double abs_frac) const{
	assert(z_ind >= 0);
	assert(z_ind < (int)grid->z.size());
	assert(dshift_l2c > 0);
	assert(lab_d >= 0);
	assert(lab_opac >= 0);
	assert(abs_frac >= 0);
	assert(abs_frac <= 1);

	// get comoving values
	double com_e = p->e * dshift_l2c;
	double com_nu = p->nu * dshift_l2c;
	double com_d  = lab_d / dshift_l2c;
	double com_opac = lab_opac / dshift_l2c;
	assert(com_e > 0);
	assert(com_nu > 0);
	assert(com_d >= 0);
	assert(com_opac >= 0);

	// set pointer to the current zone
	zone* zone;
	zone = &(grid->z[z_ind]);

	// tally in contribution to zone's distribution function (lab frame)
	// use rhat, thetahat, phihat as basis functions so rotational symmetries give accurate results
	double to_add = p->e * lab_d;
	double x=p->x[0], y=p->x[1], z=p->x[2];
	double r = sqrt(dot(p->x,p->x));
	double rp = sqrt(x*x + y*y);
	vector<double> rhat(3,0);
	vector<double> thetahat(3,0);
	vector<double> phihat(3,0);
	if(rp==0){
	    if(r==0){
	      rhat[0] = 0;
	      rhat[1] = 0;
	      rhat[2] = 1;
	    }
	    thetahat[0] = 0;
	    thetahat[1] = 1;
	    thetahat[2] = 0;
	    phihat[0] = 1;
	    phihat[1] = 0;
	    phihat[2] = 0;
	}
	else{
	  rhat[0] = x/r;
	  rhat[1] = y/r;
	  rhat[2] = z/r;
	  thetahat[0] = z/r * x/rp;
	  thetahat[1] = z/r * y/rp;
	  thetahat[2] = z/r * z/r - 1;
	  phihat[0] = -y/rp;
	  phihat[1] =  x/rp;
	  phihat[2] = 0;
	}
	vector<double> D_newbasis(3,0);
	D_newbasis[0] = dot(p->D,phihat);
	D_newbasis[1] = dot(p->D,thetahat);
	D_newbasis[2] = dot(p->D,rhat);
	normalize(D_newbasis);
	zone->distribution[p->s].count(D_newbasis, p->nu, to_add);

	// store absorbed energy in *comoving* frame (will turn into rate by dividing by dt later)
	// Extra dshift definitely needed here (two total) to convert both p->e and this_d to the comoving frame
	to_add = com_e * com_d * (com_opac*abs_frac);
	#pragma omp atomic
	zone->e_abs += to_add;
	assert(zone->e_abs >= 0);

	// store absorbed lepton number (same in both frames, except for the
	// factor of this_d which is divided out later
	double this_l_comoving = 0;
	if(species_list[p->s]->lepton_number != 0){
		this_l_comoving = com_e/(com_nu*pc::h) * com_d;
		to_add = this_l_comoving * (com_opac*abs_frac);
		if(species_list[p->s]->lepton_number == 1){
            #pragma omp atomic
			zone->nue_abs += to_add;
		}
		else if(species_list[p->s]->lepton_number == -1){
            #pragma omp atomic
			zone->anue_abs += to_add;
		}
	}
}

void transport::move(particle* p, double lab_d){
	p->x[0] += lab_d*p->D[0];
	p->x[1] += lab_d*p->D[1];
	p->x[2] += lab_d*p->D[2];
	p->t = p->t + lab_d/pc::c;
}

//--------------------------------------------------------
// Propagate a single monte carlo particle until
// it  escapes, is absorbed, or the time step ends
//--------------------------------------------------------
void transport::propagate(particle* p, const double lab_dt)
{

	ParticleEvent event;

	p->fate = moving;

	// local variables
	double lab_d = NaN;                            // distance to particle's next event
	double lab_opac = NaN, com_opac = NaN, abs_frac = NaN;                 // opacity variables
	double dshift_l2c = NaN;
	double com_nu = NaN;

	// propagate until this flag is set
	while (p->fate == moving)
	{
		assert(p->nu > 0);

		int z_ind = grid->zone_index(p->x);
		assert(z_ind >= -1);
		assert(z_ind < (int)grid->z.size());

		if(grid->good_zone(z_ind) && z_ind>=0){ // avoid handling fluff zones if unnecessary
			// doppler shift from comoving to lab (nu0/nu)
			dshift_l2c = dshift_lab_to_comoving(p,z_ind);
			assert(dshift_l2c > 0);

			// get local opacity and absorption fraction
			com_nu = p->nu * dshift_l2c;
			species_list[p->s]->get_opacity(com_nu,z_ind,&com_opac,&abs_frac);
			lab_opac = com_opac * dshift_l2c;
		}
		else{
			lab_opac = 0;
			dshift_l2c = NaN;
		}

		// decide which event happens
		which_event(p,z_ind,lab_dt,lab_opac,&lab_d,&event);
		assert(lab_d >= 0);

		// accumulate counts of radiation energy, absorption, etc
		if(grid->good_zone(z_ind) && z_ind>=0) tally_radiation(p,z_ind,dshift_l2c,lab_d,lab_opac,abs_frac);

		// move particle the distance
		move(p,lab_d);
		z_ind = grid->zone_index(p->x);

		// do the selected event
		// now the exciting bit!
		switch(event){
		    // ---------------------------------
		    // Do if interact
		    // ---------------------------------
		case interact:
			event_interact(p,z_ind,abs_frac);
			assert(p->nu > 0);
			assert(p->e > 0);
			break;

			// ---------------------------------
			// do if time step end
			// ---------------------------------
		case timeStep:
			assert(z_ind >= 0);
			p->fate = stopped;
			break;

			// ---------------------------------
			// do if crossing a boundary
			// ---------------------------------
		case boundary:{
			//it's only when the distance to move is tiny that floating point precision becomes an issue
			//keep moving by larger and larger steps until it's off grid
			{
				int i=1;
				while(z_ind>=0){
					move(p,pow(2.0,i)*lab_d);
					z_ind = grid->zone_index(p->x);
					i++;
				}
			}
			event_boundary(p,z_ind);
			break;
		}

			//-----------------------
			// nothing special happens at the zone edge
			//-----------------------
		default:
			assert(event == zoneEdge);
			assert(z_ind >= 0);
		}

		// check for core absorption
		if(p->r() < r_core) p->fate = absorbed;
	}
}
