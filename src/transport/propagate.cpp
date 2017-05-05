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

using namespace std;
namespace pc = physical_constants;

void Transport::propagate_particles()
{
	if(verbose && rank0) cout << "# Propagating particles..." << endl;

	//--- MOVE THE PARTICLES AROUND ---
	// particle list is changing size, so must go through repeatedly
	unsigned start=0, end=0;
	#pragma omp parallel
	do{
		#pragma omp single
		{
			start = end;
			end = particles.size();
		}

		#pragma omp for schedule(guided)
		for(unsigned i=start; i<end; i++){
			Particle* p = &particles[i];
			#pragma omp atomic
			n_active[p->s]++;
			if(p->fate == moving) propagate(p);
			if(p->fate == escaped){
				double nu = p->kup[3]/(2.0*pc::pi) * pc::c;
				double D[3] = {p->kup[0], p->kup[1], p->kup[2]};
				Grid::normalize_Minkowski<3>(D,3);
				#pragma omp atomic
				n_escape[p->s]++;
				#pragma omp atomic
				L_net_esc[p->s] += p->e;
				#pragma omp atomic
				N_net_esc[p->s] += p->e / (nu*pc::h);
				species_list[p->s]->spectrum.count(D,3, nu, p->e);
			}
			PRINT_ASSERT(p->fate, !=, moving);
		} //#pragma omp parallel for
	} while(particles.size()>end);

	double tot=0,core=0,fluid=0,esc=0;
	#pragma omp parallel for reduction(+:tot,core,fluid,esc)
	for(unsigned i=0; i<particles.size(); i++){
		if(particles[i].fate == moving){
			if(rank0) cout << particles[i].fate << endl;
			if(rank0) cout << i << endl;
			PRINT_ASSERT(particles[i].fate,!=,moving);
		}
		if(particles[i].fate!=rouletted) tot += particles[i].e;
		if(particles[i].fate==escaped) esc += particles[i].e;
		if(particles[i].fate==absorbed){
			if(grid->zone_index(particles[i].xup,3)>=0 ) fluid += particles[i].e;
			else core += particles[i].e;
		}
	}
	particle_total_energy += tot;
	particle_core_abs_energy += core;
	particle_fluid_abs_energy += fluid;
	particle_escape_energy += esc;

	// remove the dead particles
	particles.resize(0);
}

//--------------------------------------------------------
// Decide what happens to the particle
//--------------------------------------------------------
void Transport::which_event(LorentzHelper *lh, const int z_ind, ParticleEvent *event) const{
	PRINT_ASSERT(lh->net_opac(lab), >=, 0);
	PRINT_ASSERT(lh->p_e(lab), >, 0);
	PRINT_ASSERT(lh->p_nu(lab), >, 0);

	double d_zone     = numeric_limits<double>::infinity();
	double d_interact = numeric_limits<double>::infinity();
	double d_boundary = numeric_limits<double>::infinity();
	PRINT_ASSERT(z_ind, >=, -1);

	if(z_ind >= 0){ //i.e. within the simulation region
		// FIND D_ZONE= ====================================================================
		// maximum step size inside zone
		double d_zone_min = step_size * grid->zone_min_length(z_ind);
		double d_zone_boundary = (1.0+TINY) * grid->zone_cell_dist(lh->p_xup(), z_ind);
		d_zone = max(d_zone_min, d_zone_boundary);
		PRINT_ASSERT(d_zone, >, 0);

		// FIND D_INTERACT =================================================================
		double relevant_opacity = lh->tau_opac(lab);
		if (relevant_opacity == 0) d_interact = numeric_limits<double>::infinity();
		else               d_interact = static_cast<double>(lh->p_tau() / relevant_opacity);
		PRINT_ASSERT(d_interact,>=,0);
	}

	// FIND D_BOUNDARY ================================================================
	d_boundary = grid->lab_dist_to_boundary(lh);
	PRINT_ASSERT(d_boundary, >, 0);

	// find out what event happens (shortest distance) =====================================
	*event  = zoneEdge;
	double d_smallest = d_zone;
	if( d_interact <= d_smallest ){
		*event = interact;
		d_smallest = d_interact;
	}
	if( d_boundary <= d_smallest ){
		*event = boundary;
		d_smallest = d_boundary;
		if(z_ind >= 0) d_smallest *= (1.0 + TINY); // bump just over the boundary if in simulation domain
		else           d_smallest *= (1.0 - TINY); // don't overshoot outward through the inner boundary
	}
	lh->set_distance<lab>(d_smallest);
	PRINT_ASSERT(lh->distance(lab), >=, 0);
	PRINT_ASSERT(lh->distance(lab), <, INFINITY);
}


void Transport::event_boundary(LorentzHelper *lh, const int z_ind) const{
	PRINT_ASSERT(z_ind==-1, ||, z_ind==-2);

	// if outside the domain
	if(z_ind == -2){
		int new_ind = z_ind;
		if(reflect_outer){
			grid->reflect_outer(lh);
			PRINT_ASSERT(lh->p_fate(), ==, moving);
			new_ind = grid->zone_index(lh->p_xup(),3);
			PRINT_ASSERT(new_ind, >=, 0);
			PRINT_ASSERT(new_ind, <, (int)grid->z.size());
			PRINT_ASSERT(lh->p_nu(lab), >, 0);
		}
		else{
			grid->symmetry_boundaries(lh);
			new_ind = grid->zone_index(lh->p_xup(),3);
		}

		if(new_ind < 0) lh->set_p_fate(escaped);
	}

	// if inside the inner boundary
	if(z_ind==-1){
		if(grid->radius(lh->p_xup(),3) < r_core) lh->set_p_fate(absorbed);
		else PRINT_ASSERT(lh->p_fate(), ==, moving); // the particle just went into the inner boundary
	}
}

void Transport::distribution_function_basis(const double D[3], const double xyz[3], double D_newbasis[3]) const{
	double x=xyz[0], y=xyz[1], z=xyz[2];
	double r = sqrt(Grid::dot_Minkowski<3>(xyz,xyz,3));
	double rp = sqrt(x*x + y*y);
	double rhat[3]     = {0,0,0};
	double thetahat[3] = {0,0,0};
	double phihat[3]   = {0,0,0};
	if(rp==0){
		rhat[2] = z>0 ? -1.0 : 1.0;
	    thetahat[1] = 1;
	    phihat[0] = 1;
	}
	else{
		rhat[0] = x/r;
		rhat[1] = y/r;
		rhat[2] = z/r;
		thetahat[0] = z/r * x/rp;
		thetahat[1] = z/r * y/rp;
		thetahat[2] = -rp / r;
		phihat[0] = -y/rp;
		phihat[1] =  x/rp;
		phihat[2] = 0;
	}
	D_newbasis[0] = Grid::dot_Minkowski<3>(D,thetahat,3);
	D_newbasis[1] = Grid::dot_Minkowski<3>(D,phihat,3);
	D_newbasis[2] = Grid::dot_Minkowski<3>(D,rhat,3);
}
void Transport::tally_radiation(const LorentzHelper *lh, const int z_ind) const{
	PRINT_ASSERT(z_ind, >=, 0);
	PRINT_ASSERT(z_ind, <, (int)grid->z.size());
	PRINT_ASSERT(lh->distance(lab), >=, 0);
	PRINT_ASSERT(lh->abs_fraction(), >=, 0);
	PRINT_ASSERT(lh->abs_fraction(), <=, 1);
	PRINT_ASSERT(lh->p_e(com), >, 0);
	PRINT_ASSERT(lh->p_nu(com), >, 0);
	PRINT_ASSERT(lh->distance(com), >=, 0);
	PRINT_ASSERT(lh->net_opac(com), >=, 0);
	double to_add = 0;

	// set pointer to the current zone
	Zone* zone;
	zone = &(grid->z[z_ind]);

	// get the distribution function basis to use
	// use rhat, thetahat, phihat as basis functions so rotational symmetries give accurate results (if distribution_polar_basis)
	double Dlab[3], D_newbasis[3];
	lh->p_D(lab,Dlab,3);
	if(distribution_polar_basis) distribution_function_basis(Dlab,lh->p_xup(),D_newbasis);
	else for(int i=0;i<3;i++) D_newbasis[i] = Dlab[i];
	Grid::normalize_Minkowski<3>(D_newbasis,3);

	// tally in contribution to zone's distribution function (lab frame)
	if(exponential_decay) to_add = lh->p_e(lab) / lh->abs_opac(lab) * (1.0 - exp(-lh->abs_opac(lab) * lh->distance(lab)));
	else to_add = lh->p_e(lab) * lh->distance(lab);
	PRINT_ASSERT(to_add,<,INFINITY);
	zone->distribution[lh->p_s()]->count(D_newbasis, 3, lh->p_nu(com), to_add);

	// store absorbed energy in *comoving* frame (will turn into rate by dividing by dt later)
	if(exponential_decay) to_add = lh->p_e(com) * (1.0 - exp(-lh->abs_opac(com) * lh->distance(com)));
	else to_add = lh->p_e(com) * lh->distance(com) * (lh->abs_opac(com));
	#pragma omp atomic
	zone->e_abs += to_add;
	PRINT_ASSERT(zone->e_abs, >=, 0);

	// store absorbed lepton number (same in both frames, except for the
	// factor of this_d which is divided out later
	if(species_list[lh->p_s()]->lepton_number != 0){
		to_add /= (lh->p_nu(com)*pc::h);
		if(species_list[lh->p_s()]->lepton_number == 1){
            #pragma omp atomic
			zone->nue_abs += to_add;
		}
		else if(species_list[lh->p_s()]->lepton_number == -1){
            #pragma omp atomic
			zone->anue_abs += to_add;
		}
	}
}

void Transport::move(LorentzHelper *lh, int *z_ind){
	PRINT_ASSERT(lh->p_tau(),>=,0);
	PRINT_ASSERT(lh->distance(lab),>=,0);
	PRINT_ASSERT(lh->net_opac(lab),>=,0);

	// translate the particle
	grid->integrate_geodesic(lh);

	// reduce the particle's remaining optical depth
	if(lh->tau_opac(lab)>0){
		double old_tau = lh->p_tau();
		double new_tau = lh->p_tau() - lh->tau_opac(lab) * lh->distance(lab);
		PRINT_ASSERT(new_tau,>=,-TINY*old_tau);
		lh->set_p_tau( max(0.0,new_tau) );
	}

	// appropriately reduce the particle's energy
	if(exponential_decay){
		lh->scale_p_e( exp(-lh->abs_opac(lab) * lh->distance(lab)) );
		window(lh,*z_ind);
	}

	// re-evaluate the particle's zone index
	*z_ind = grid->zone_index(lh->p_xup(),3);

}

void Transport::get_opacity(LorentzHelper *lh, const int z_ind) const{
	if(grid->good_zone(z_ind) && z_ind>=0){ // avoid handling fluff zones if unnecessary

		// get local opacity and absorption fraction
		double a=-1,s=-1;
		species_list[lh->p_s()]->get_opacity(lh->p_nu(com), z_ind, &a, &s);
		lh->set_opac<com>(a,s);
	}
	else{
		lh->set_opac<com>(0,0);
	}
	PRINT_ASSERT(lh->net_opac(lab),>=,0);
	PRINT_ASSERT(lh->net_opac(lab),<,INFINITY);
}
//--------------------------------------------------------
// Propagate a single monte carlo particle until
// it  escapes, is absorbed, or the time step ends
//--------------------------------------------------------
void Transport::propagate(Particle* p)
{

	ParticleEvent event;
	LorentzHelper lh(exponential_decay);
	lh.set_p<lab>(p);

	PRINT_ASSERT(lh.p_fate(), ==, moving);

	while (lh.p_fate() == moving)
	{
		PRINT_ASSERT(lh.p_nu(lab), >, 0);

		int z_ind = grid->zone_index(lh.p_xup(),3);
		PRINT_ASSERT(z_ind, >=, -1);
		PRINT_ASSERT(z_ind, <, (int)grid->z.size());

		// set up the LorentzHelper
		double v[3];
		grid->interpolate_fluid_velocity(lh.p_xup(),3,v,3,z_ind);
		lh.set_v(v,3);

		// get all the opacities
		get_opacity(&lh,z_ind);

		// decide which event happens
		which_event(&lh,z_ind,&event);

		// accumulate counts of radiation energy, absorption, etc
		if(grid->good_zone(z_ind) && z_ind>=0) tally_radiation(&lh,z_ind);

		// move particle the distance
		move(&lh, &z_ind);
		if(event != boundary) PRINT_ASSERT(lh.p_tau(), >=, -TINY*(lh.distance(lab) * lh.tau_opac(lab)));

		// do the selected event
		// now the exciting bit!
		if(lh.p_fate() != rouletted) switch(event){
		    // ---------------------------------
		    // Do if interacting with the fluid
		    // ---------------------------------
		case interact:
			PRINT_ASSERT(lh.distance(lab) * lh.net_opac(lab), >=, lh.p_tau());
			event_interact(&lh,z_ind);
			if(lh.p_fate()==moving){
				PRINT_ASSERT(lh.p_nu(lab), >, 0);
				PRINT_ASSERT(lh.p_e(lab), >, 0);
			}

			break;

			// ---------------------------------
			// do if crossing a boundary
			// ---------------------------------
		case boundary:{
			//it's only when the distance to move is tiny that floating point precision becomes an issue
			//keep moving by larger and larger steps until it's off grid
			{
				z_ind = grid->zone_index(lh.p_xup(),3);
				int i=1;
				while(z_ind>=0){
					double tweak_distance = TINY*grid->zone_min_length(z_ind) * i*i;
					PRINT_ASSERT(tweak_distance,<,2*grid->zone_min_length(z_ind));
					lh.set_p_tau(lh.p_tau() + tweak_distance*lh.tau_opac(lab) ); // a hack to prevent tau<0
					lh.set_distance<lab>(tweak_distance);
					move(&lh, &z_ind);

					i++;
				}
			}
			event_boundary(&lh,z_ind);
			break;
		}

			//-----------------------
			// nothing special happens at the zone edge
			//-----------------------
		default:
			PRINT_ASSERT(event, ==, zoneEdge);
			PRINT_ASSERT(z_ind, >=, 0);
		}

		// check for core absorption
		if(grid->radius(lh.p_xup(),3) < r_core) lh.set_p_fate(absorbed);
	}

	// copy particle back out of LorentzHelper
	*p = lh.particle_copy(lab);
	PRINT_ASSERT(p->fate, !=, moving);
}
