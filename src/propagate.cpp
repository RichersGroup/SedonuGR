#pragma warning disable 161
#include <limits>
#include <omp.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <cassert>
#include "transport.h"
#include "physical_constants.h"
#include "species_general.h"

namespace pc = physical_constants;

void transport::propagate_particles(const double dt)
{
	vector<int> n_active(species_list.size(),0);
	vector<int> n_escape(species_list.size(),0);
	double e_esc = 0;
    #pragma omp parallel shared(n_active,n_escape,e_esc) firstprivate(dt)
	{

		//--- MOVE THE PARTICLES AROUND ---
        #pragma omp for schedule(guided) reduction(+:e_esc)
		for(int i=0; i<particles.size(); i++){
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
			if(steady_state && verbose) cout << "# e_esc = " << e_esc << "erg" << endl;
			if(dt<0) assert(particles.size()==0);
		}

	} //#pragma omp parallel

	//--- OUPUT ESCAPE STATISTICS ---
	if (verbose){
		for(int i=0; i<species_list.size(); i++){
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
	double d_zone     = numeric_limits<double>::infinity();
	double d_time     = numeric_limits<double>::infinity();
	double d_interact = numeric_limits<double>::infinity();
	double d_boundary = numeric_limits<double>::infinity();
	double tau_r;                               // random optical depth
	double opac_lab;                            // opacity in the lab frame
	assert(p->ind >= -1);
	assert(opac >=0);

	if(p->ind >= 0){ //i.e. within the simulation region
		// set pointer to current zone
		zone* zone;
		zone = &(grid->z[p->ind]);

		// FIND D_ZONE= ====================================================================
		// maximum step size inside zone
		d_zone = step_size * grid->zone_min_length(p->ind);

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
	if(dt>0){
		double tstop = t_now + dt;
		d_time = (tstop - p->t)*pc::c;
	}

	// FIND D_BOUNDARY ================================================================
	d_boundary = grid->dist_to_boundary(p);

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
		if(p->ind >= 0) *d_smallest *= (1.0 + grid->tiny); // bump just over the boundary if in simulation domain
		else *d_smallest *= (1.0 - grid->tiny);            // don't overshoot outward through the inner boundary
	}
}


//--------------------------------------------------------
// Propagate a single monte carlo particle until
// it  escapes, is absorbed, or the time step ends
//--------------------------------------------------------
void transport::propagate(particle* p, const double dt) const
{
	assert(p->ind >= -1);
	assert(p->ind < grid->z.size());
	ParticleEvent event;

	p->fate = moving;

	// local variables
	double this_d = 0;                            // distance to particle's next event
	double dshift = 0;                                   // doppler shift
	double opac = 0, abs_frac = 0;                 // opacity variables
	double this_E=0, this_E_comoving=0, this_l_comoving=0; // for calculating radiation energy and energy/lepton number absorbed
	double z=0,z2=0;                                     // random numbers

	// propagate until this flag is set
	while (p->fate == moving)
	{
		assert(p->ind >= -1);
		// set pointer to current zone
		zone* zone;
		zone = &(grid->z[p->ind]);

		// doppler shift from comoving to lab
		dshift = dshift_comoving_to_lab(p);
		// get local opacity and absorption fraction
		species_list[p->s]->get_opacity(p,dshift,&opac,&abs_frac);

		// decide which event happens
		which_event(p,dt,dshift,opac,&this_d,&event);

		// tally in contribution to zone's radiation energy (both *lab* frame)
		this_E = p->e*this_d;
        #pragma omp atomic
		zone->e_rad += this_E;

		// store absorbed energy in *comoving* frame (will turn into rate by dividing by dt later)
		// Extra dshift definitely needed here (two total) to convert both p->e and this_d to the comoving frame
		this_E_comoving = this_E * dshift * dshift;
        #pragma omp atomic
		zone->e_abs += this_E_comoving * (opac*abs_frac*zone->eps_imc);

		// store absorbed lepton number (same in both frames, except for the
		// factor of this_d which is divided out later
		if(species_list[p->s]->lepton_number != 0){
			this_l_comoving = species_list[p->s]->lepton_number * p->e/(p->nu*pc::h) * this_d*dshift;
            #pragma omp atomic
			zone->l_abs += this_l_comoving * (opac*abs_frac*zone->eps_imc);
		}

		// TODO - put back in radiation force tally here
		// fx_rad =

		// move particle the distance
		p->x[0] += this_d*p->D[0];
		p->x[1] += this_d*p->D[1];
		p->x[2] += this_d*p->D[2];

		// advance the time
		p->t = p->t + this_d/pc::c;

		// get zone location now
		p->ind = grid->get_zone(p->x);

		// now the exciting bit!
		switch(event){
		// ---------------------------------
		// Do if interact
		// ---------------------------------
		case interact:
			assert(p->ind >= 0);
			// random number to check for scattering or absorption
			z = rangen.uniform();

			// decide whether to scatter or absorb
			if (z > abs_frac) isotropic_scatter(p,0); // actual scatter - do not redistribute energy
			else
			{
				// if this is an iterative calculation, radiative equilibrium is always assumed.
				if(radiative_eq) isotropic_scatter(p,1);          // particle lives, energy redistributed
				else{
					z2 = rangen.uniform();
					if (z2 > zone->eps_imc) isotropic_scatter(p,1); // particle lives, energy redistributed
					else p->fate = absorbed;                        // particle dies
				}
			}
			break;

		// ---------------------------------
		// do if time step end
		// ---------------------------------
		case timeStep:
			assert(p->ind >= 0);
			p->fate = stopped;
			break;

		// ---------------------------------
		// do if crossing a boundary
		// ---------------------------------
		case boundary:
			if(p->ind >= 0){
				cout << p->ind << endl;
				cout << p->r()-1.0e7 << endl;
				cout << p->x_dot_d() << endl;
				cout << grid->get_zone(p->x) << endl;
			}
			assert(p->ind==-1 || p->ind==-2);

			// if outside the domain
			if(p->ind == -2){
				if(reflect_outer){
					grid->reflect_outer(p);
					assert(p->fate == moving);
					assert(p->ind >= 0);
				}
				else p->fate = escaped;
			}

			// if inside the inner boundary
			if(p->ind==-1){
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
			if(p->ind < 0){
				cout << p->r() << endl;
				cout << this_d << endl;
			}
			assert(event == zoneEdge);
			assert(p->ind >= 0);
		}

		// check for core absorption
		if(p->r() < r_core) p->fate = absorbed;
	}
}
