#include <omp.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include "transport.h"
#include "species_general.h"
#include "global_options.h"

void transport::propagate_particles(const double lab_dt)
{
	if(verbose && rank0) cout << "# Propagating particles..." << endl;
	double e_esc_lab = 0;

	//--- MOVE THE PARTICLES AROUND ---
    #pragma omp parallel for schedule(guided) reduction(+:e_esc_lab)
	for(unsigned i=0; i<particles.size(); i++){
		particle* p = &particles[i];
        #pragma omp atomic
		n_active[p->s]++;
		propagate(p,lab_dt);
		if(p->fate == escaped){
			#pragma omp atomic
			n_escape[p->s]++;
			species_list[p->s]->spectrum.count(p->t, p->nu, p->e, p->D);
			e_esc_lab += p->e;
		}
	} //#pragma omp parallel fpr

	// store the escaped energy in the transport class
	L_esc_lab += e_esc_lab; // to be normalized later

	// remove the dead particles
	remove_dead_particles();
}

void transport::remove_dead_particles(){
	vector<particle>::iterator pIter = particles.begin();
	while(pIter != particles.end()){
		if(pIter->fate==absorbed || pIter->fate==escaped){
			*pIter = particles[particles.size()-1];
			particles.pop_back();
		}
		else pIter++;
	}
	if(steady_state) assert(particles.empty());
}


//--------------------------------------------------------
// Decide what happens to the particle
//--------------------------------------------------------
void transport::which_event(const particle *p, const double dt, const double lab_opac,
		double *d_smallest, ParticleEvent *event) const{
	assert(lab_opac >= 0);

	const int z_ind = grid->zone_index(p->x);
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
	if(!steady_state){
		assert(dt > 0);
		double tstop = t_now + dt;
		d_time = (tstop - p->t)*pc::c;
		assert(d_time > 0);
	}

	// FIND D_BOUNDARY ================================================================
	d_boundary = grid->lab_dist_to_boundary(p);
	assert(d_boundary >= 0);

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


void transport::do_event(const ParticleEvent event, const double abs_frac, particle* p){
	// get zone location now
	int z_ind = grid->zone_index(p->x);
	assert(z_ind >= -2);
	assert(z_ind < (int)grid->z.size());

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
	case boundary:
		event_boundary(p,z_ind);
		break;

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

void transport::event_boundary(particle* p, const int z_ind) const{
	assert(z_ind==-1 || z_ind==-2);

	// if outside the domain
	if(z_ind == -2){
		if(reflect_outer){
			grid->reflect_outer(p);
			assert(p->fate == moving);
			assert(grid->zone_index(p->x) >= 0);
			assert(grid->zone_index(p->x) < (int)grid->z.size());
			assert(p->nu > 0);
		}
		else p->fate = escaped;
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

void transport::tally_radiation(const particle* p, const double dshift_l2c, const double lab_d, const double lab_opac, const double abs_frac) const{
	int z_ind = grid->zone_index(p->x);
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
	double to_add=0;

	// tally in contribution to zone's radiation energy (both *lab* frame)
	to_add = com_e * com_d;
	#pragma omp atomic
	zone->e_rad += to_add;
	assert(zone->e_rad >= 0);

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
		this_l_comoving = species_list[p->s]->lepton_number * com_e/(com_nu*pc::h) * com_d;
		to_add = this_l_comoving * (com_opac*abs_frac);
        #pragma omp atomic
		zone->l_abs += to_add;
	}

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
		int on_grid = (z_ind >= 0);

		if(grid->good_zone(z_ind) && on_grid){ // avoid handling fluff zones if unnecessary
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
		which_event(p,lab_dt,lab_opac,&lab_d,&event);
		assert(lab_d >= 0);

		// accumulate counts of radiation energy, absorption, etc
		if(grid->good_zone(z_ind) && on_grid) tally_radiation(p,dshift_l2c,lab_d,lab_opac,abs_frac);

		// move particle the distance
		p->x[0] += lab_d*p->D[0];
		p->x[1] += lab_d*p->D[1];
		p->x[2] += lab_d*p->D[2];
		p->t = p->t + lab_d/pc::c;

		// do the selected event
		do_event(event,abs_frac,p);
	}
}
