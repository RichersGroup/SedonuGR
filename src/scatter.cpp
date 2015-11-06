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
#include <math.h>
#include <gsl/gsl_rng.h>
#include "species_general.h"
#include "transport.h"

//------------------------------------------------------------
// physics of absorption/scattering
//------------------------------------------------------------
void transport::event_interact(particle* p, const int z_ind, const double abs_frac, const double lab_opac){
	assert(z_ind >= 0);
	assert(z_ind < (int)grid->z.size());
	assert(abs_frac >= 0.0);
	assert(abs_frac <= 1.0);
	assert(p->e > 0);

	int do_absorb_partial  = !radiative_eq;
	int do_absorb_reemit   =  radiative_eq;

	// particle is transformed to the comoving frame
	transform_lab_to_comoving(p,z_ind);

	// absorb part of the packet
	if(do_absorb_partial){
		if(abs_frac < 1.0){
			p->e *= (1.0 - abs_frac);
			isotropic_scatter(p);
		}
		if(p->e<min_packet_energy) window(p,z_ind);
		assert(p->e > 0.0);
	}

	// absorb the particle and let the fluid re-emit another particle
	else if(do_absorb_reemit){
		re_emit(p,z_ind);
		L_net_lab[p->s] += p->e;
		assert(p->e > 0.0);
	}

	// particle is transformed back to the lab frame
	assert(p->e > 0);
	transform_comoving_to_lab(p,z_ind);

	// resample the path length
	sample_tau(p,z_ind,lab_opac);

	// sanity checks
	assert(p->nu > 0);
	assert(p->e > 0);
}


// decide whether to kill a particle
void transport::window(particle* p, const int z_ind){
	// Roulette if too low energy
	while(p->e<min_packet_energy && p->fate!=absorbed){
		if(rangen.uniform() < 0.5) p->fate = absorbed;
		else p->e *= 2.0;
	}
	// split if too high energy, if enoug space, and if in important region
	while(p->e>max_packet_energy && particles.size()<max_particles && species_list[p->s]->importance(p->nu,z_ind)>=1.0){
		p->e /= 2.0;
		particle pnew = *p;
		#pragma omp critical
		particles.push_back(pnew);
	}

	if(!particles.size()<max_particles && verbose && rank0)
		cout << "WARNING: max_particles is too small to allow splitting." << endl;
}

// re-emission, done in COMOVING frame
void transport::re_emit(particle* p, const int z_ind) const{
	assert(z_ind >= 0);
	assert(z_ind < (int)grid->z.size());

	// reset the particle properties
	isotropic_scatter(p);
	p->s = sample_zone_species(z_ind);
	species_list[p->s]->sample_zone_nu(*p,z_ind);

	// tally into zone's emitted energy
	grid->z[z_ind].e_emit += p->e;

	// sanity checks
	assert(p->nu > 0);
	assert(p->s >= 0);
	assert(p->s < (int)species_list.size());
}

// isotropic scatter, done in COMOVING frame
void transport::isotropic_scatter(particle* p) const
{
	// Randomly generate new direction isotropically in comoving frame
	double mu  = 1 - 2.0*rangen.uniform();
	double phi = 2.0*pc::pi*rangen.uniform();
	double smu = sqrt(1 - mu*mu);
	p->D[0] = smu*cos(phi);
	p->D[1] = smu*sin(phi);
	p->D[2] = mu;
	normalize(p->D);
}

//---------------------------------------------------------------------
// Randomly select an optical depth through which a particle will move.
//---------------------------------------------------------------------
void transport::sample_tau(particle *p, const int z_ind, const double lab_opac) const{
	assert(opt_depth_bias>=0);
	double taubar = lab_opac * grid->zone_min_length(z_ind);
	double tau_r = opt_depth_bias * taubar * -1.0*log(1.0 - rangen.uniform());
	p->e *= exp(tau_r * (1./(opt_depth_bias*taubar) - 1.0));
}



