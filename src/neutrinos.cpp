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
  std::vector<double>stg = lua->vector<double>("nut_spec_time_grid");
  std::vector<double>sng = lua->vector<double>("nut_spec_log_nu_grid");
  int nmu  = lua->scalar<int>("nut_n_mu");
  int nphi = lua->scalar<int>("nut_n_phi");
  spectrum.init(stg,sng,nmu,nphi);
  spectrum.set_name("neutrino_spectrum.dat");

  // read opacity parameters
  grey_opac = lua->scalar<double>("nut_grey_opacity");
  eps       = lua->scalar<double>("nut_epsilon");
  if(nulibID == 0){
    lepton_number = 1;
    name = "Electron Neutrinos";}
  else if(nulibID == 1){
    lepton_number = -1;
    name = "Electron Anti-Neutrinos";}
  else if(nulibID == 2){
    lepton_number = 0;
    if(num_nut_species == 3)      name = "Mu/Tau Anti/Neutrinos";
    else if(num_nut_species == 4) name = "Mu/Tau Neutrinos";
    else if(num_nut_species == 6) name = "Mu Neutrinos";
    else name = "ERROR";}
  else if(nulibID == 3){
    lepton_number = 0;
    if(num_nut_species == 4)      name = "Mu/Tau Anti-Neutrinos";
    else if(num_nut_species == 6) name = "Mu Antineutrino";
    else name = "ERROR";}
  else if(nulibID == 4){
    lepton_number = 0;
    if(num_nut_species == 6) name = "Tau Neutrinos";
    else name = "ERROR";}
  else if(nulibID == 5){
    lepton_number = 0;
    if(num_nut_species == 6) name = "Tau Anti-Neutrinos";
    else name = "ERROR";}
  else{
    cout << "ERROR: Sedona does not know how to deal with a neutrino ID of " << nulibID << "." << endl;
    exit(16);}

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

  // set neutrino's min and max values
  T_min  =  nulib_get_Tmin();
  T_max  =  nulib_get_Tmax();
  Ye_min =  nulib_get_Yemin();
  Ye_max =  nulib_get_Yemax();
}


//-----------------------------------------------------------------
// set emissivity, abs. opacity, and scat. opacity in zones
//-----------------------------------------------------------------
void neutrinos::set_eas(int zone_index)
{
    zone* z = &(sim->grid->z[zone_index]);

    z->eps_imc = 1;
    nulib_get_eas_arrays(z->rho, z->T_gas, z->Ye, nulibID,
			 emis[zone_index], abs_opac[zone_index], scat_opac[zone_index]);
    emis[zone_index].normalize();
}
