#pragma warning disable 161
#include <vector>
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
  std::vector<double>stg = lua->vector<double>("spec_time_grid");
  std::vector<double>sng = lua->vector<double>("spec_nu_grid");
  int nmu  = lua->scalar<int>("spec_n_mu");
  int nphi = lua->scalar<int>("spec_n_phi");
  spectrum.init(stg,sng,nmu,nphi);
  spectrum.set_name("neutrino_spectrum.dat");

  // read opacity parameters
  grey_opac = lua->scalar<double>("nut_grey_opacity");
  eps       = lua->scalar<double>("nut_epsilon");

  // set lepton number
  if(num_nut_species==3){
	  if(nulibID == 0)   lepton_number =  1;
	  if(nulibID == 1)   lepton_number = -1;
	  if(nulibID == 2)   lepton_number =  0;
  }
  else{
	  if(nulibID%2 == 0) lepton_number =  1;
	  if(nulibID%2 == 1) lepton_number = -1;
  }

  // set names and spectrum weight
  double weight;
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

  // read in the frequency table
  nulib_get_nu_grid(nu_grid);

  // allocate space for the grid eas spectrum containers
  abs_opac.resize(sim->grid->z.size());
  scat_opac.resize(sim->grid->z.size());
  emis.resize(sim->grid->z.size());

  // now allocate space for each eas spectrum
  if(sim->do_core) core_emis.resize(nu_grid.size());
  for (int i=0; i<emis.size();      i++)      emis[i].resize(nu_grid.size());
  for (int i=0; i<abs_opac.size();  i++)  abs_opac[i].resize(nu_grid.size());
  for (int i=0; i<scat_opac.size(); i++) scat_opac[i].resize(nu_grid.size());

  // set up core neutrino emission spectrum function (erg/s)
  // normalize to core luminosity. constants don't matter.
  if(sim->do_core){
    double T_core = lua->scalar<double>("T_core");
    double L_core = lua->scalar<double>("L_core");
    double chem_pot = 0;
    #pragma omp parallel for ordered
    for (int j=0;j<nu_grid.size();j++)
    {
      double nu  = nu_grid.center(j);
      double dnu = nu_grid.delta(j);
      double bb  = nu*nu*nu*fermi_dirac(T_core,chem_pot,nu)*dnu;
      #pragma omp ordered
      core_emis.set_value(j,bb);
    }
    core_emis.normalize();
    core_emis.N = weight * L_core / 6.0;
  }

  // set neutrino's min and max values
  T_min  =  nulib_get_Tmin();
  T_max  =  nulib_get_Tmax();
  Ye_min =  nulib_get_Yemin();
  Ye_max =  nulib_get_Yemax();
  rho_min = nulib_get_rhomin();
  rho_max = nulib_get_rhomax();
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

//-----------------------------------------------------------------
// Calculate the fermi-dirac function (erg/s/cm^2/Hz/ster)
// (normalized to 1, not total luminosity)
//-----------------------------------------------------------------
double neutrinos::fermi_dirac(const double T, const double chem_pot, const double nu) const
{
	double zeta = (pc::h*nu - chem_pot)/pc::k/T;
	return 1.0 / (exp(zeta) + 1);
}
