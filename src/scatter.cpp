#include <math.h>
#include <gsl/gsl_rng.h>
#include "species_general.h"
#include "transport.h"
#include "global_options.h"

//------------------------------------------------------------
// physics of isotropic scattering
//------------------------------------------------------------
void transport::isotropic_scatter(particle* p, const int redistribute) const
{
	int z_ind = grid->zone_index(p->x);
	assert(z_ind >= 0);
	assert(z_ind < (int)grid->z.size());

	// transform into the comoving frame
	transform_lab_to_comoving(p);

	// Randomly generate new direction isotropically in comoving frame
	double mu  = 1 - 2.0*rangen.uniform();
	double phi = 2.0*pc::pi*rangen.uniform();
	double smu = sqrt(1 - mu*mu);
	p->D[0] = smu*cos(phi);
	p->D[1] = smu*sin(phi);
	p->D[2] = mu;
	p->normalize_direction();

	// change wavelength and species, if needed
	if (redistribute){
		p->s = sample_zone_species(z_ind);
		p->nu = species_list[p->s]->sample_zone_nu(z_ind);
		assert(p->nu > 0);
		assert(p->s >= 0);
		assert(p->s < species_list.size());
	}

	// lorentz transform back to lab frame
	transform_comoving_to_lab(p);

	// sanity checks
	assert(p->nu > 0);
}

