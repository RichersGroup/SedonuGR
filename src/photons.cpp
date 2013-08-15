#include "photons.h"
#include "transport.h"
#include <vector>
#include "physical_constants.h"

namespace pc = physical_constants;
using namespace std;

//----------------------------------------------------------------
// called from species_general::init (photon-specific stuff)
//----------------------------------------------------------------
void photons::myInit(Lua* lua)
{
  // intialize output spectrum
  std::vector<double>stg = lua->vector<double>("spec_time_grid");
  std::vector<double>sng = lua->vector<double>("spec_nu_grid");
  int nmu  = lua->scalar<int>("n_mu");
  int nphi = lua->scalar<int>("n_phi");
  spectrum.init(stg,sng,nmu,nphi);
  spectrum.set_name("optical_spectrum.dat");

  // read opacity parameters
  gray_abs_opac  = lua->scalar<double>("gray_abs_opacity");
  gray_scat_opac = lua->scalar<double>("gray_scat_opacity");
  eps            = lua->scalar<double>("epsilon");
  double nu_start = lua->scalar<double>("nu_start");
  double nu_stop  = lua->scalar<double>("nu_stop");
  int      n_nu   = lua->scalar<int>("n_nu");

  // initialize the  frequency grid
  nu_grid.init(nu_start,nu_stop,n_nu);

  // allocate space for the grid eas spectrum containers
  if(gray_abs_opac <=0)  abs_opac.resize(sim->grid->z.size());
  if(gray_scat_opac<=0) scat_opac.resize(sim->grid->z.size());
  emis.resize(sim->grid->z.size());

  // now allocate space for each eas spectrum
  core_emis.resize(nu_grid.size());
  for (int i=0; i<abs_opac.size();  i++)  abs_opac[i].resize(nu_grid.size());
  for (int i=0; i<scat_opac.size(); i++) scat_opac[i].resize(nu_grid.size());
  for (int i=0; i<emis.size();      i++)      emis[i].resize(nu_grid.size());

  // set up core emission spectrum function (now a blackbody) 
  double T_core = lua->scalar<double>("T_core");
  for (int j=0;j<nu_grid.size();j++)
  {
    double nu  = nu_grid.center(j);
    double dnu = nu_grid.delta(j);
    double bb  = blackbody_nu(T_core,nu)*dnu;
    core_emis.set_value(j,bb); 
  }
  core_emis.normalize();
}
