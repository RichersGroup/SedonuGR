#include <omp.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include "transport.h"
#include "species_general.h"
#include "global_options.h"

void transport::propagate_particles(const double dt)
{
	vector<int> n_active(species_list.size(),0);
	vector<int> n_escape(species_list.size(),0);
	double e_esc = 0;
#pragma omp parallel shared(n_active,n_escape,e_esc) firstprivate(dt)
	{

		//--- MOVE THE PARTICLES AROUND ---
#pragma omp for schedule(guided) reduction(+:e_esc)
		for(unsigned i=0; i<particles.size(); i++){
			particle* p = &particles[i];
#pragma omp atomic
			n_active[p->s]++;
			propagate(p,dt);
			if(p->fate == escaped){
#pragma omp atomic
				n_escape[p->s]++;
				species_list[p->s]->spectrum.count(p->t, p->nu, p->e, p->D);
				e_esc += p->e;
			}
		} //implied barrier

		//--- REMOVE THE DEAD PARTICLES ---
#pragma omp single
		{
			vector<particle>::iterator pIter = particles.begin();
			while(pIter != particles.end()){
				if(pIter->fate==absorbed || pIter->fate==escaped){
					*pIter = particles[particles.size()-1];
					particles.pop_back();
				}
				else pIter++;
			}
		}

		// report energy escaped
#pragma omp single
		{
			if(rank0 && verbose) cout << "# e_esc = " << e_esc << " erg" << endl;
			if(steady_state) assert(particles.size()==0);
		}

	} //#pragma omp parallel

	//--- OUPUT ESCAPE STATISTICS ---
	if (rank0 && verbose){
		for(unsigned i=0; i<species_list.size(); i++){
			double per_esc = (100.0*n_escape[i])/n_active[i];
			if(n_active[i]>0) printf("# %i/%i %s escaped. (%3.2f%%)\n", n_escape[i], n_active[i], species_list[i]->name.c_str(), per_esc);
			else printf("# No active %s.\n", species_list[i]->name.c_str());
		}
	}
}


//--------------------------------------------------------
// Decide what happens to the particle
//--------------------------------------------------------
void transport::which_event(const particle *p, const double dt, const double dshift, const double opac,
		double *d_smallest, ParticleEvent *event) const{
	assert(dshift > 0);
	assert(opac >= 0);

	const int z_ind = grid->zone_index(p->x);
	double d_zone     = numeric_limits<double>::infinity();
	double d_time     = numeric_limits<double>::infinity();
	double d_interact = numeric_limits<double>::infinity();
	double d_boundary = numeric_limits<double>::infinity();
	double tau_r = NaN;                               // random optical depth
	double opac_lab = NaN;                            // opacity in the lab frame
	assert(z_ind >= -1);

	if(z_ind >= 0){ //i.e. within the simulation region
		// FIND D_ZONE= ====================================================================
		// maximum step size inside zone
		d_zone = step_size * grid->zone_min_length(z_ind);
		assert(d_zone > 0);

		// FIND D_INTERACT =================================================================
		// convert opacity from comoving to lab frame for the purposes of
		// determining the interaction distance in the lab frame
		// This corresponds to equation 90.8 in Mihalas&Mihalas. You multiply
		// the comoving opacity by nu_0 over nu, which is why you
		// multiply by dshift instead of dividing by dshift here
		opac_lab = opac*dshift;
		// random optical depth to next interaction
		tau_r = -1.0*log(1.0 - rangen.uniform());
		// step size to next interaction event
		d_interact  = tau_r/opac_lab;
		if (opac_lab == 0) d_interact = numeric_limits<double>::infinity();
		assert(d_interact>0);
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
	d_boundary = grid->dist_to_boundary(p);
	assert(d_boundary >= 0);

	// find out what event happens (shortest distance)
	*d_smallest = numeric_limits<double>::infinity();
	if( d_interact <= *d_smallest ){
		*event  = interact;
		*d_smallest = d_interact;
	}
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
}


//--------------------------------------------------------
// Propagate a single monte carlo particle until
// it  escapes, is absorbed, or the time step ends
//--------------------------------------------------------
void transport::propagate(particle* p, const double dt) const
{

	ParticleEvent event;

	p->fate = moving;

	// local variables
	double this_d = NaN;                            // distance to particle's next event
	double opac = NaN, abs_frac = NaN;                 // opacity variables
	double rand = NaN;

	// propagate until this flag is set
	while (p->fate == moving)
	{
		int z_ind = grid->zone_index(p->x);
		assert(z_ind >= -1);
		assert(z_ind < (int)grid->z.size());

		assert(p->nu > 0);
		// set pointer to current zone
		zone* zone;
		zone = &(grid->z[z_ind]);

		// doppler shift from comoving to lab
		double dshift = dshift_comoving_to_lab(p);
		assert(dshift > 0);

		// get local opacity and absorption fraction
		species_list[p->s]->get_opacity(p,z_ind,dshift,&opac,&abs_frac);

		// decide which event happens
		which_event(p,dt,dshift,opac,&this_d,&event);
		assert(this_d >= 0);

		// tally in contribution to zone's radiation energy (both *lab* frame)
		double this_E = p->e*this_d;
		assert(this_E > 0);
#pragma omp atomic
		zone->e_rad += this_E;

		// store absorbed energy in *comoving* frame (will turn into rate by dividing by dt later)
		// Extra dshift definitely needed here (two total) to convert both p->e and this_d to the comoving frame
		double this_E_comoving = this_E * dshift * dshift;
#pragma omp atomic
		zone->e_abs += this_E_comoving * (opac*abs_frac);

		// store absorbed lepton number (same in both frames, except for the
		// factor of this_d which is divided out later
		double this_l_comoving = 0;
		if(species_list[p->s]->lepton_number != 0){
			this_l_comoving = species_list[p->s]->lepton_number * p->e/(p->nu*pc::h) * this_d*dshift;
#pragma omp atomic
			zone->l_abs += this_l_comoving * (opac*abs_frac);
		}

		// move particle the distance
		p->x[0] += this_d*p->D[0];
		p->x[1] += this_d*p->D[1];
		p->x[2] += this_d*p->D[2];

		// advance the time
		p->t = p->t + this_d/pc::c;

		// get zone location now
		z_ind = grid->zone_index(p->x);
		assert(z_ind < (int)grid->z.size());

		// now the exciting bit!
		switch(event){
		// ---------------------------------
		// Do if interact
		// ---------------------------------
		case interact:
			assert(z_ind >= 0);
			// random number to check for scattering or absorption
			rand = rangen.uniform();

			// decide whether to scatter or absorb
			if (rand > abs_frac) isotropic_scatter(p,0); // actual scatter - do not redistribute energy
			else{
				// if this is an iterative calculation, radiative equilibrium is always assumed.
				if(radiative_eq) isotropic_scatter(p,1);          // particle lives, energy redistributed
				else p->fate = absorbed;
			}
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
				if(p->x_dot_d() >= 0){
					// set the particle just outside the inner boundary
					cout << "ERROR: have not yet implemented passing out through the inner boundary without overshooting" << endl;
					exit(5);
				}
				else assert(p->fate == moving);
			}
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
}
