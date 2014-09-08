#include <math.h>
#include <gsl/gsl_rng.h>
#include "species_general.h"
#include "transport.h"
#include "global_options.h"

//------------------------------------------------------------
// physics of absorption/scattering
//------------------------------------------------------------
void transport::event_interact(particle* p, const int z_ind, const double abs_frac){
	assert(z_ind >= 0);
	assert(z_ind < (int)grid->z.size());
	const double rand = rangen.uniform();
	int do_absorb_kill     = (rand<abs_frac) && (!radiative_eq);
	int do_absorb_reemit   = (rand<abs_frac) && ( radiative_eq);
	int do_elastic_scatter = (rand>abs_frac);

	// transform into the comoving frame
	if(do_absorb_kill) p->fate = absorbed;
	else{
		transform_lab_to_comoving(p);
		if(do_elastic_scatter) isotropic_scatter(p);
		if(do_absorb_reemit)   re_emit(p,z_ind);
		transform_comoving_to_lab(p);

		// tally the re-emitted energy in the lab frame
		if(do_absorb_reemit) L_net_lab += p->e;

		// sanity checks
		assert(p->nu > 0);
		assert(p->e > 0);
	}

}

// re-emission, done in COMOVING frame
void transport::re_emit(particle* p, const int z_ind) const{
	assert(z_ind >= 0);
	assert(z_ind < (int)grid->z.size());

	// reset the particle properties
	isotropic_scatter(p);
	p->s = sample_zone_species(z_ind);
	p->nu = species_list[p->s]->sample_zone_nu(z_ind);

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
	p->normalize_direction();
}



