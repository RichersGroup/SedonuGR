#include "neutrinos.h"
#include "transport.h"
#include <vector>
#include "physical_constants.h"
#include "nulib_interface.h"

namespace pc = physical_constants;
using namespace std;

//----------------------------------------------------------------
// called from species_general::init (photon-specific stuff)
//----------------------------------------------------------------
void neutrinos::myInit(Lua* lua)
{
  // intialize output spectrum
  std::vector<double>stg = lua->vector<double>("spec_time_grid");
  std::vector<double>sng = lua->vector<double>("spec_nu_grid");
  int nmu  = lua->scalar<int>("n_mu");
  int nphi = lua->scalar<int>("n_phi");
  spectrum.init(stg,sng,nmu,nphi);
  spectrum.set_name("neutrino_spectrum.dat");

  // read opacity parameters
  gray_abs_opac  = lua->scalar<double>("gray_abs_opacity");
  gray_scat_opac = lua->scalar<double>("gray_scat_opacity");
  eps            = lua->scalar<double>("epsilon");
  if(nulibID == 0) electron_number = 1;
  else if (nulibID == 1) electron_number = -1;
  else electron_number = 0;

  // read in the frequency table
  nulib_get_nu_grid(nu_grid);

  // allocate space for the grid eas spectrum containers
  abs_opac.resize(sim->grid->z.size());
  scat_opac.resize(sim->grid->z.size());
  emis.resize(sim->grid->z.size());

  // now allocate space for each eas spectrum
  core_emis.resize(nu_grid.size());
  for (int i=0; i<emis.size();      i++)      emis[i].resize(nu_grid.size());
  for (int i=0; i<abs_opac.size();  i++)  abs_opac[i].resize(nu_grid.size());
  for (int i=0; i<scat_opac.size(); i++) scat_opac[i].resize(nu_grid.size());

  // set up core neutrino emission spectrum function
  double rho_core = lua->scalar<double>("rho_core");
  double T_core   = lua->scalar<double>("T_core");
  double Ye_core  = lua->scalar<double>("Ye_core");
  //just place holders so we can use the function to get core_emis
  vector<double> tmp1(nu_grid.size(),0), tmp2(nu_grid.size(),0);
  nulib_get_eas_arrays(rho_core, T_core, Ye_core, nulibID,
		       core_emis, tmp1, tmp2);
  core_emis.normalize();
}


//-----------------------------------------------------------------
// set emissivity, abs. opacity, and scat. opacity in zones
//-----------------------------------------------------------------
void neutrinos::set_eas()
{
  for (int i=0; i < sim->grid->z.size(); i++)
  {
    zone* z = &(sim->grid->z[i]);
    nulib_get_eas_arrays(z->rho, z->T_gas, z->Ye, nulibID,
			 emis[i], abs_opac[i], scat_opac[i]);
    emis[i].normalize();
  }
}
