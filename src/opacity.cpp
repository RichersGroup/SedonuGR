#include <gsl/gsl_rng.h>
#include <math.h>
#include "transport.h"
#include "physical_constants.h"
#include "Lua.h"
#include "photons.h"

namespace pc = physical_constants;

//-----------------------------------------------------------------
// get opacity at the frequency
//-----------------------------------------------------------------
void species_general::get_eas(particle &p, double dshift, double* e, double* a, double* s)
{
  // comoving frame frequency
  double nu = p.nu*dshift;

  // pointer to current zone
  zone *zone = &(sim->grid->z[p.ind]);

  // TODO - better way to output variables consistently? (e is not the emissivity)
  // emissivity
  if(eps >= 0) *e = eps;
  else *e = nu_grid.value_at<cdf_array>(nu,emis[p.ind]);

  // absorption opacity
  if(gray_abs_opac >= 0) *a = zone->rho * gray_abs_opac;
  else *a = nu_grid.value_at< vector<double> >(nu,abs_opac[p.ind]);

  // scattering opacity
  if(gray_scat_opac >= 0) *s = zone->rho * gray_scat_opac;
  else *s = nu_grid.value_at< vector<double> >(nu,scat_opac[p.ind]);
}
