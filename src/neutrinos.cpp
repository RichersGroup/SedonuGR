#pragma warning disable 161
#include <vector>
#include <cassert>
#include <limits>
#include "neutrinos.h"
#include "transport.h"
#include "Lua.h"
#include "physical_constants.h"
#include "nulib_interface.h"

namespace pc = physical_constants;
using namespace std;

//----------------------------------------------------------------
// called from species_general::init (neutrino-specific stuff)
//----------------------------------------------------------------
void neutrinos::myInit(Lua* lua)
{
  // intialize output spectrum
  std::vector<double>stg = lua->vector<double>("nut_spec_time_grid");
  std::vector<double>sng = lua->vector<double>("nut_spec_nu_grid");
  int nmu  = lua->scalar<int>("nut_spec_n_mu");
  int nphi = lua->scalar<int>("nut_spec_n_phi");
  spectrum.init(stg,sng,nmu,nphi);
  spectrum.set_name("neutrino_spectrum.dat");

  // set lepton number
  if(nulibID == 0)   lepton_number =  1;
  if(nulibID == 1)   lepton_number = -1;
  if(nulibID >= 2)   lepton_number =  0;

  // set names and spectrum weight
  if     (nulibID == 0) {name = "Electron Neutrinos";              weight = 1.0;}
  else if(nulibID == 1) {name = "Electron Anti-Neutrinos";         weight = 1.0;}
  else if(nulibID == 2){
    if     (num_nut_species == 3) {name = "Mu/Tau Anti/Neutrinos"; weight = 4.0;}
    else if(num_nut_species == 4) {name = "Mu/Tau Neutrinos";      weight = 2.0;}
    else if(num_nut_species == 6) {name = "Mu Neutrinos";          weight = 1.0;}
    else                          {name = "ERROR";                 weight = -1;}}
  else if(nulibID == 3){
    if     (num_nut_species == 4) {name = "Mu/Tau Anti-Neutrinos"; weight = 2.0;}
    else if(num_nut_species == 6) {name = "Mu Antineutrino";       weight = 1.0;}
    else                          {name = "ERROR";                 weight = -1;}}
  else if(nulibID == 4){
    if(num_nut_species == 6)      {name = "Tau Neutrinos";         weight = 1.0;}
    else                          {name = "ERROR";                 weight = -1;}}
  else if(nulibID == 5){
    if(num_nut_species == 6)      {name = "Tau Anti-Neutrinos";    weight = 1.0;}
    else                          {name = "ERROR";                 weight = -1;}}
  else{
    cout << "ERROR: Sedona does not know how to deal with a neutrino ID of " << nulibID << "." << endl;
    exit(16);}

  // set up the frequency table
  grey_opac     = lua->scalar<double>("nut_grey_opacity");
  grey_abs_frac = lua->scalar<double>("nut_grey_abs_frac");
  if(grey_opac < 0){
	  nulib_get_nu_grid(nu_grid);

	  // set neutrino's min and max values
	  T_min  =  nulib_get_Tmin();
	  T_max  =  nulib_get_Tmax();
	  Ye_min =  nulib_get_Yemin();
	  Ye_max =  nulib_get_Yemax();
	  rho_min = nulib_get_rhomin();
	  rho_max = nulib_get_rhomax();

  }
  else{
	  double nu_start = lua->scalar<double>("nut_nugrid_start");
	  double nu_stop  = lua->scalar<double>("nut_nugrid_stop");
	  int      n_nu   = lua->scalar<int>("nut_nugrid_n");
	  nu_grid.init(nu_start,nu_stop,n_nu);

	  // set neutrino's min and max values
	  T_min   =  1.0;
	  T_max   =  1e12;
	  Ye_min  = 0;
	  Ye_max  = 1;
	  rho_min = -numeric_limits<double>::infinity();
	  rho_max =  numeric_limits<double>::infinity();

  }

}


//-----------------------------------------------------------------
// set emissivity, abs. opacity, and scat. opacity in zones
//-----------------------------------------------------------------
void neutrinos::set_eas(int zone_index)
{
	zone* z = &(sim->grid->z[zone_index]);

	if(grey_opac < 0){ // get opacities and emissivity from NuLib
		nulib_get_eas_arrays(z->rho, z->T_gas, z->Ye, nulibID,
				emis[zone_index], abs_opac[zone_index], scat_opac[zone_index]);
		emis[zone_index].normalize();
	}

	else{ // get emissivity from blackbody and the grey opacity
		assert(grey_abs_frac>=0 && grey_abs_frac<=1.0);
		for (int j=0;j<nu_grid.size();j++)
		{
			double nu  = nu_grid.center(j);        // (Hz)
			double dnu = nu_grid.delta(j);         // (Hz)
			double bb  = blackbody(z->T_gas,0*pc::MeV_to_ergs,nu)*dnu;  // (erg/s/cm^2/ster)

			double a = grey_opac*z->rho*grey_abs_frac;
			double s = grey_opac*z->rho*(1.0-grey_abs_frac);

			emis[zone_index].set_value(j,a*bb); // (erg/s/cm^3/ster)
			abs_opac[zone_index][j] = a;        // (1/cm)
			scat_opac[zone_index][j] = s;       // (1/cm)
		}
		emis[zone_index].normalize();
	}
}


//-----------------------------------------------------------------
// Calculate the fermi-dirac blackbody function (erg/s/cm^2/Hz/ster)
//-----------------------------------------------------------------
double neutrinos::blackbody(const double T, const double chem_pot, const double nu) const
{
  double zeta = (pc::h*nu - chem_pot)/pc::k/T;
  return (pc::h/pc::c/pc::c) * nu*nu*nu / (exp(zeta) + 1.0);
}
