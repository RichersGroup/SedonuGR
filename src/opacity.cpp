#include <gsl/gsl_rng.h>
#include <math.h>
#include "transport.h"
#include "physical_constants.h"
#include "Lua.h"

namespace pc = physical_constants;

//-----------------------------------------------------------------
// allocate room for the opacities we will store
//-----------------------------------------------------------------
void transport::initialize_opacity(Lua *lua)
{
  // read opacity parameters
  this->grey_opac = lua->scalar<double>("grey_opacity");
  double nu_start = lua->scalar<double>("nu_start");
  double nu_stop  = lua->scalar<double>("nu_stop");
  int      n_nu   = lua->scalar<int>("n_nu");

  // initialize the  frequency grid
  nu_grid.init(nu_start,nu_stop,n_nu);
  //  if (verbose) nu_grid.print();

  // allocate space for the opac/emis vectors
  for (int i=0;i<grid->n_zones;i++)
  {
    grid->z[i].opac.resize(nu_grid.size());
    grid->z[i].emis.resize(nu_grid.size());
  }

  // allocate space for any core emission spectrum
  this->n_inject  = lua->scalar<int>("n_inject");
  this->r_core    = lua->scalar<double>("r_core");
  this->core_emis.resize(nu_grid.size());
}


//-----------------------------------------------------------------
// set opacity and emissivities
//-----------------------------------------------------------------
void transport::set_opacity()
{
  for (int i=0;i<grid->n_zones;i++)
  {
    // fleck factors
    //double Tg    = grid->z[i].T_gas;
    //double fleck_beta=4.0*pc::a*pow(Tg,4)/(grid->z[i].e_gas*grid->z[i].rho);
    //double tfac  = pc::c*grey_opac*epsilon*grid->z[i].rho*t_step;
    //double f_imc = fleck_alpha*fleck_beta*tfac;
    //grid->z[i].eps_imc = 1.0/(1.0 + f_imc);
    //if (radiative_eq) grid->z[i].eps_imc = 1.;

    grid->z[i].eps_imc = 1;

    zone* z = &(grid->z[i]);
    for (int j=0;j<nu_grid.size();j++)
    {
      double nu  = nu_grid.value(j);
      double bb  = blackbody_nu(z->T_gas,nu);
      z->opac[j] = grey_opac*z->rho;
      z->emis.set_value(j,grey_opac*bb);
    }
    z->emis.normalize();
  }
}



//-----------------------------------------------------------------
// get opacity at the frequency
//-----------------------------------------------------------------
void transport::get_opacity(particle &p, double dshift, double &opac, double &eps)
{
  // comoving frame frequency
  double nu = p.nu*dshift;

  // pointer to current zone
  zone *zone = &(grid->z[p.ind]);

  // get opacity if it is an optical photon. 
  if (p.type == photon)
  {
    // interpolate opacity at the local comving frame frequency
    opac = nu_grid.value_at(nu,zone->opac);
    eps  = this->epsilon;
    
    // grey opacity flag will now override, if set
    if (this->grey_opac > 0)
      { opac = zone->rho*this->grey_opac; eps = this->epsilon; }
  }
}


//-----------------------------------------------------------------
// Klein_Nishina correction to the Compton cross-section
// assumes energy x is in MeV
//-----------------------------------------------------------------
double transport::klein_nishina(double x)
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
double transport::blackbody_nu(double T, double nu)
{
  double zeta = pc::h*nu/pc::k/T;
  return 2.0*nu*nu*nu*pc::h/pc::c/pc::c/(exp(zeta)-1);
}
