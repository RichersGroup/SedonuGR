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
void species_general::get_opacity(particle &p, double dshift, double* opac, double* abs_frac)
{
  // comoving frame frequency
  double nu = p.nu*dshift;

  // absorption opacity
  double a = 0;
  if(grey_opac<0) a = nu_grid.value_at< vector<double> >(nu,abs_opac[p.ind]);
  else a = sim->grid->z[p.ind].rho * grey_opac;

  // scattering opacity
  double s = 0;
  if(eps<0 && grey_opac<0) s = nu_grid.value_at< vector<double> >(nu,scat_opac[p.ind]);

  // output - net opacity
  *opac = a+s;

  // output - absorption fraction
  if(eps < 0) *abs_frac = a/(a+s);
  else        *abs_frac = eps;
}
