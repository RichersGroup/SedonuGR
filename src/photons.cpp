#pragma warning disable 161
#include <limits>
#include <vector>
#include "photons.h"
#include "transport.h"
#include "Lua.h"
#include "physical_constants.h"

namespace pc = physical_constants;
using namespace std;

//----------------------------------------------------------------
// called from species_general::init (photon-specific stuff)
//----------------------------------------------------------------
void photons::myInit(Lua* lua)
{
  // set name
  name = "Photons";

  // poison unused zone properties
  #pragma omp parallel for
  for(int i=0; i<sim->grid->z.size(); i++) sim->grid->z[i].Ye = -1.0e99;

  // intialize output spectrum
  std::vector<double>stg = lua->vector<double>("spec_time_grid");
  std::vector<double>sng = lua->vector<double>("spec_nu_grid");
  int nmu  = lua->scalar<int>("spec_n_mu");
  int nphi = lua->scalar<int>("spec_n_phi");
  spectrum.init(stg,sng,nmu,nphi);
  spectrum.set_name("optical_spectrum.dat");

  // read opacity parameters
  grey_opac       = lua->scalar<double>("grey_opacity");
  eps             = lua->scalar<double>("epsilon");
  double nu_start = lua->scalar<double>("nu_start");
  double nu_stop  = lua->scalar<double>("nu_stop");
  int      n_nu   = lua->scalar<int>("n_nu");
  lepton_number   = 0;

  // initialize the  frequency grid
  nu_grid.init(nu_start,nu_stop,n_nu);

  // allocate space for the grid eas spectrum containers
  abs_opac.resize(sim->grid->z.size());
  scat_opac.resize(sim->grid->z.size());
  emis.resize(sim->grid->z.size());

  // now allocate space for each eas spectrum
  if(sim->n_emit_core > 0) core_emis.resize(nu_grid.size());
  #pragma omp parallel for
  for (int i=0; i<abs_opac.size();  i++){
    abs_opac[i].resize(nu_grid.size());
    scat_opac[i].resize(nu_grid.size());
    emis[i].resize(nu_grid.size());
  }

  // set up core emission spectrum function (now a blackbody) (erg/s)
  // normalized to core luminosity. constants don't matter.
  if(sim->n_emit_core > 0){
    double T_core = lua->scalar<double>("T_core");
    double L_core = lua->scalar<double>("L_core");
    #pragma omp parallel for ordered
    for (int j=0;j<nu_grid.size();j++)
    {
      double nu  = nu_grid.center(j);
      double dnu = nu_grid.delta(j);
      double bb  = planck(T_core,nu)*dnu;
      #pragma omp ordered
      core_emis.set_value(j,bb); 
    }
    core_emis.normalize();
    core_emis.N = L_core;
  }

  // set photon's min and max values
  T_min   =  1.0;
  T_max   =  1e12;
  Ye_min  = -numeric_limits<double>::infinity();
  Ye_max  =  numeric_limits<double>::infinity();
  rho_min = -numeric_limits<double>::infinity();
  rho_max =  numeric_limits<double>::infinity();
}


//-----------------------------------------------------------------
// set emissivity, abs. opacity, and scat. opacity in zones
//-----------------------------------------------------------------
void photons::set_eas(int zone_index)
{
  zone* z = &(sim->grid->z[zone_index]);
  if(grey_opac >= 0)
  {
    // fleck factors
    //double Tg    = sim->grid->z[zone_index].T_gas;
    //double fleck_beta=4.0*pc::a*pow(Tg,4)/(sim->grid->z[zone_index].e_gas*sim->grid->z[zone_index].rho);
    //double tfac  = pc::c*grey_opac*epsilon*sim->grid->z[zone_index].rho*t_step;
    //double f_imc = fleck_alpha*fleck_beta*tfac;
    //sim->grid->z[zone_index].eps_imc = 1.0/(1.0 + f_imc);
    //if (sim->radiative_eq) sim->grid->z[zone_index].eps_imc = 1.;
    z->eps_imc = 1;

    // leave serial. Parrallelized threads call this function.
    for (int j=0;j<nu_grid.size();j++)
    {
      double nu  = nu_grid.center(j);        // (Hz)
      double dnu = nu_grid.delta(j);         // (Hz)
      double bb  = planck(z->T_gas,nu)*dnu;  // (erg/s/cm^2/ster)
      emis[zone_index].set_value(j,grey_opac*eps*bb*z->rho); // (erg/s/cm^3/Hz/ster)
    }
    emis[zone_index].normalize();
  }
}



//-----------------------------------------------------------------
// Klein_Nishina correction to the Compton cross-section
// assumes energy x is in MeV
//-----------------------------------------------------------------
double photons::klein_nishina(const double x_input) const
{
  // divide by m_e c^2 = 0.511 MeV
  double x = x_input/pc::m_e_MeV;
  double logfac = log(1 + 2*x);
  double term1 = (1+x)/x/x/x*(2*x*(1+x)/(1+2*x) - logfac);
  double term2 = 1.0/2.0/x*logfac;
  double term3 = -1.0*(1 + 3*x)/(1+2*x)/(1+2*x);
  double KN    = .75*(term1 + term2 + term3);
  return KN;
}

// calculate planck function (erg/s/cm^2/Hz/ster)
double photons::planck(const double T, const double nu) const
{
  double zeta = pc::h*nu/pc::k/T;
  return 2.0*nu*nu*nu*pc::h/pc::c/pc::c/(exp(zeta)-1);
}
