#include <gsl/gsl_rng.h>
#include <math.h>
#include "transport.h"
#include "physical_constants.h"
#include "Lua.h"
#include "photons.h"

namespace pc = physical_constants;

//-----------------------------------------------------------------
// set opacity and emissivities
//-----------------------------------------------------------------
void photons::set_eas()
{
  for (int i=0;i<sim->grid->z.size();i++)
  {
    // fleck factors
    //double Tg    = sim->grid->z[i].T_gas;
    //double fleck_beta=4.0*pc::a*pow(Tg,4)/(sim->grid->z[i].e_gas*sim->grid->z[i].rho);
    //double tfac  = pc::c*grey_opac*epsilon*sim->grid->z[i].rho*t_step;
    //double f_imc = fleck_alpha*fleck_beta*tfac;
    //sim->grid->z[i].eps_imc = 1.0/(1.0 + f_imc);
    //if (sim->radiative_eq) sim->grid->z[i].eps_imc = 1.;
    sim->grid->z[i].eps_imc = 1;

    zone* z = &(sim->grid->z[i]);
    for (int j=0;j<nu_grid.size();j++)
    {
      double nu  = nu_grid.x[j];
      double bb  = blackbody_nu(z->T_gas,nu);
      abs_opac[i][j] = gray_abs_opac*z->rho;
      emis[i].set_value(j,gray_abs_opac*bb);
    }
    emis[i].normalize();
  }
}



//-----------------------------------------------------------------
// get opacity at the frequency
//-----------------------------------------------------------------
void species_general::get_eas(particle &p, double dshift, double &e, double &a, double &s)
{
  // comoving frame frequency
  double nu = p.nu*dshift;

  // pointer to current zone
  zone *zone = &(sim->grid->z[p.ind]);

  // TODO - check units
  // emissivity
  if(eps > 0) e = eps;
  else e = nu_grid.value_at<cdf_array>(nu,emis[p.ind]);

  // absorption opacity
  if(gray_abs_opac > 0) a = gray_abs_opac;
  else a = nu_grid.value_at< vector<double> >(nu,abs_opac[p.ind]);
    
  // scattering opacity
  if(gray_scat_opac > 0) s = gray_scat_opac;
  else s = nu_grid.value_at< vector<double> >(nu,scat_opac[p.ind]);
}


//-----------------------------------------------------------------
// Klein_Nishina correction to the Compton cross-section
// assumes energy x is in MeV
//-----------------------------------------------------------------
double photons::klein_nishina(double x)
{
  // divide by m_e c^2 = 0.511 MeV
  x = x/pc::m_e_MeV;
  double logfac = log(1 + 2*x);
  double term1 = (1+x)/x/x/x*(2*x*(1+x)/(1+2*x) - logfac);
  double term2 = 1.0/2.0/x*logfac;
  double term3 = -1.0*(1 + 3*x)/(1+2*x)/(1+2*x);
  double KN    = .75*(term1 + term2 + term3);
  return KN;
}

// calculate planck function in frequency units
double photons::blackbody_nu(double T, double nu)
{
  double zeta = pc::h*nu/pc::k/T;
  return 2.0*nu*nu*nu*pc::h/pc::c/pc::c/(exp(zeta)-1);
}
