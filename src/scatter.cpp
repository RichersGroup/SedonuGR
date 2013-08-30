#include "species_general.h"
#include <math.h>
#include <gsl/gsl_rng.h>
#include "transport.h"
#include "physical_constants.h"

namespace pc = physical_constants;



//------------------------------------------------------------
// physics of isotropic scattering
//------------------------------------------------------------
void transport::isotropic_scatter(particle &p, int re_emit)
{

  // get doppler shift from lab to comoving frame
  double dshift_in = dshift_lab_to_comoving(p);
  
  // transform energy, frequency into comoving frame
  p.e  *= dshift_in;
  p.nu *= dshift_in;

  // Randomly generate new direction isotropically in comoving frame
  double mu  = 1 - 2.0*gsl_rng_uniform(rangen);
  double phi = 2.0*pc::pi*gsl_rng_uniform(rangen);
  double smu = sqrt(1 - mu*mu);
  p.D[0] = smu*cos(phi);
  p.D[1] = smu*sin(phi);
  p.D[2] = mu;
  
  // change wavelength and species, if needed
//  int s;
//  if (re_emit){
//	  s = sim->sample_zone_species(p.ind);
//	  p.nu = sim->species_list[s]->sample_zone_nu(p.ind);
//  }

  // lorentz transform back to lab frame
  transform_comoving_to_lab(p);

  //
}
