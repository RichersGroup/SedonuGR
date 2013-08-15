#include "species_general.h"
#include <math.h>
#include <gsl/gsl_rng.h>
#include "transport.h"
#include "physical_constants.h"

namespace pc = physical_constants;



//------------------------------------------------------------
// physics of isotropic scattering
//------------------------------------------------------------
void species_general::isotropic_scatter(particle &p, int redistribute)
{

  // get doppler shift from lab to comoving frame
  double dshift_in = dshift_lab_to_comoving(p);
  
  // transform energy, frequency into comoving frame
  p.e  *= dshift_in;
  p.nu *= dshift_in;

  // Randomly generate new direction isotropically in comoving frame
  double mu  = 1 - 2.0*gsl_rng_uniform(sim->rangen);
  double phi = 2.0*pc::pi*gsl_rng_uniform(sim->rangen);
  double smu = sqrt(1 - mu*mu);
  p.D[0] = smu*cos(phi);
  p.D[1] = smu*sin(phi);
  p.D[2] = mu;
  
  // change wavelength, if needed
  if (redistribute) p.nu = sample_zone_nu(p.ind);

  // lorentz transform back to lab frame
  transform_comoving_to_lab(p);
}
