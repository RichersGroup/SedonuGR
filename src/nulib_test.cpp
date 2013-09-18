#include <iomanip>
#include <vector>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string.h>
#include "zone.h"
#include "nulib_interface.h"
#include "physical_constants.h"

namespace pc = physical_constants;
using namespace std;

int main(int argc, char* argv[]){
  using namespace std;
  if(argc!=7){
    cout << "Usage: nulib_test.exe path_to_nulib_table.h5 rho(g/cm^3) T(MeV) Ye E(MeV) nulibID" << endl;
    exit(1);
  }

  //set test inputs
  real rho        = atof(argv[2]);           // g/cm^3
  real T          = atof(argv[3])/pc::k_MeV; // K
  real ye         = atof(argv[4]);
  double myenergy = atof(argv[5]);           // MeV
  double myfreq   = myenergy     /pc::h_MeV; // Hz
  int    nulibID  = atoi(argv[6]);
  
  //read in the nulib table
  cout << "initializing nulib" << endl;
  string filename = argv[1];
  nulib_init(filename);
  
  //=============//
  // GRID OUTPUT //
  //=============//
  ofstream grid_file;
  grid_file.open("grids.dat");

  locate_array nu_grid; // Hz
  nulib_get_nu_grid(nu_grid);
  grid_file << "Energy Grid (MeV)" << endl;
  for(int i=0; i<nu_grid.size(); i++) grid_file << nu_grid.x[i]*pc::h_MeV << endl;
  grid_file << endl;

  vector<real> ye_grid;
  nulib_get_Ye_array(ye_grid);
  grid_file << "Ye grid" << endl;
  for(int i=0; i<ye_grid.size(); i++) grid_file << ye_grid[i] << endl;
  grid_file << endl;
  
  vector<real> T_grid; // K
  nulib_get_T_array(T_grid);
  grid_file << "T grid (MeV)" << endl;
  for(int i=0; i<T_grid.size(); i++) grid_file << T_grid[i]*pc::k_MeV << endl;
  grid_file << endl;

  vector<real> rho_grid; // g/cm^3
  nulib_get_rho_array(rho_grid);
  grid_file << "rho grid (g/cm^3)" << endl;
  for(int i=0; i<rho_grid.size(); i++) grid_file << rho_grid[i] << endl;
  grid_file << endl;

  grid_file.close();

  // read in the number of species and groups in the table
  int ns = nulib_get_nspecies();
  int nrho = rho_grid.size();
  int nye = ye_grid.size();
  int nT = T_grid.size();
  int ng = nu_grid.size();
  cout << ns << " species and " << ng << " groups" << endl;

  //make vectors of appropriate sizes
  cout << "making vectors" << endl;
  vector<double> absopac  (ng,0); // cm^-1
  vector<double> scatopac (ng,0); // cm^-1
  vector<double> pure_emis(ng,0); // erg/cm^3/s/ster/Hz
  cdf_array emis;                 // erg/cm^3/s/ster
  emis.resize(ng);

  //===================//
  // SINGLE LINE PLOTS //
  //===================//
  double e,a,s;

  ofstream eas_rho;
  eas_rho.open("eas_rho.dat");
  cout << "generating eas_rho.dat" << endl;
  eas_rho << "# T(MeV):" << T*pc::k_MeV << " ye:" << ye << " E(MeV):" << myenergy << endl; 
  eas_rho << setw(25) << "# rho(g/cm^3)" << setw(25) << "emis(int)(erg/cm^3/s/ster)" << setw(25) << "absopac(cm^2/g)"<< setw(25) <<"scatopac(cm^2/g)" << endl;
  for(int j=0; j<nrho; j++){
    nulib_get_eas_arrays(rho_grid[j], T, ye, nulibID, emis, absopac, scatopac);
    e = emis.get(nu_grid.size()-1);
    a = nu_grid.value_at(myfreq, absopac);
    s = nu_grid.value_at(myfreq, scatopac);
    eas_rho << setw(25) << rho_grid[j] << setw(25) << e << setw(25) << a/rho_grid[j] << setw(25) << s/rho_grid[j] << endl;
  }

  ofstream eas_T;
  eas_T.open("eas_T.dat");
  cout << "generating eas_T.dat" << endl;
  eas_T << "# rho(g/cm^3):" << rho << " ye:" << ye << " E(MeV):" << myenergy << endl; 
  eas_T << setw(25) << "# T(MeV)" << setw(25) << "emis(int)(erg/cm^3/s/ster)" << setw(25) << "absopac(cm^2/g)" <<setw(25) << "scatopac(cm^2/g)" << endl;
  for(int j=0; j<nT; j++){
    nulib_get_eas_arrays(rho, T_grid[j], ye, nulibID, emis, absopac, scatopac);
    e = emis.get(nu_grid.size()-1);
    a = nu_grid.value_at(myfreq, absopac);
    s = nu_grid.value_at(myfreq, scatopac);
    eas_T <<setw(25) << T_grid[j]*pc::k_MeV <<setw(25) << e <<setw(25) << a/rho <<setw(25) << s/rho << endl;
  }

  ofstream eas_ye;
  eas_ye.open("eas_ye.dat");
  cout << "generating eas_ye.dat" << endl;
  eas_ye << "# rho(g/cm^3):" << rho << " T(MeV):" << T*pc::k_MeV << " E(MeV):" << myenergy << endl; 
  eas_ye << setw(25) << "# ye" << setw(25) << "emis(int)(erg/cm^3/s/ster)" << setw(25) << "absopac(cm^2/g)" << setw(25) << "scatopac(cm^2/g)" << endl;
  for(int j=0; j<nye; j++){
    nulib_get_eas_arrays(rho, T, ye_grid[j], nulibID, emis, absopac, scatopac);
    e = emis.get(nu_grid.size()-1);
    a = nu_grid.value_at(myfreq, absopac);
    s = nu_grid.value_at(myfreq, scatopac);
    eas_ye <<setw(25) << ye_grid[j] <<setw(25) << e <<setw(25) << a/rho <<setw(25) << s/rho << endl;
  }

  eas_rho.close();
  eas_T.close();
  eas_ye.close();

  //=====================//
  // MULTIPLE LINE PLOTS //
  //=====================//
  // variation with rho
  cout << "generating eas_E_rho.dat" << endl;
  ofstream eas_E_rho ("eas_E_rho.dat" );

  eas_E_rho << "# T(MeV)="<<T*pc::k_MeV << " ye="<<ye << " rho_min(g/cm^3)="<<nulib_get_rhomin() << " rho_max(g/cm^3)="<<nulib_get_rhomax() << " Nrho=" << rho_grid.size() << endl;
  eas_E_rho << setw(25) << "# rho(g/cm^3)" << setw(25) << "E(MeV)" << setw(25) << "emis(erg/cm^3/s/ster/Hz)" << setw(25) << "absopac(cm^2/g)" << setw(25) << "scatopac(cm^2/g)" << endl;

  for(int i=0; i<rho_grid.size(); i++){
    nulib_get_eas_arrays(rho_grid[i], T, ye, nulibID, emis, absopac, scatopac);
    nulib_get_pure_emis (rho_grid[i], T, ye, nulibID, pure_emis);
    for(int j=0; j<nu_grid.size(); j++){
      e = pure_emis[j];
      a = nu_grid.value_at(nu_grid[j], absopac);
      s = nu_grid.value_at(nu_grid[j], scatopac);
      eas_E_rho << setw(25) << rho_grid[i] <<setw(25) << nu_grid[j]*pc::h_MeV <<setw(25) << e  <<setw(25) << a/rho_grid[j] <<setw(25) << s/rho_grid[j] << endl;
    }
    eas_E_rho << endl;
  }

  eas_E_rho.close();


  // variation with T
  cout << "generating eas_E_T.dat" << endl;
  ofstream eas_E_T ("eas_E_T.dat" );

  eas_E_T << "# rho(g/cm^3)="<<rho << " ye="<<ye << " T_min(MeV)="<<nulib_get_Tmin()*pc::k_MeV << " T_max(MeV)="<<nulib_get_Tmax()*pc::k_MeV << " NT=" << T_grid.size() << endl;
  eas_E_T << setw(25) << "# T(MeV)" << setw(25) << "E(MeV)" << setw(25) << "emis(erg/cm^3/s/ster/Hz)" << setw(25) << "absopac(cm^2/g)" << setw(25) << "scatopac(cm^2/g)" << endl;

  for(int i=0; i<T_grid.size(); i++){
    nulib_get_eas_arrays(rho, T_grid[i], ye, nulibID, emis, absopac, scatopac);
    nulib_get_pure_emis (rho_grid[i], T, ye, nulibID, pure_emis);
    for(int j=0; j<nu_grid.size(); j++){
      e = pure_emis[j];
      a = nu_grid.value_at(nu_grid[j], absopac);
      s = nu_grid.value_at(nu_grid[j], scatopac);
      eas_E_T << setw(25) << T_grid[i]*pc::k_MeV << setw(25) << nu_grid[j]*pc::h_MeV << setw(25) << e << setw(25) << a/rho << setw(25) << s/rho << endl;
    }
    eas_E_T << endl;
  }

  eas_E_T.close();


  // variation with ye
  cout << "generating eas_E_ye.dat" << endl;
  ofstream eas_E_ye ("eas_E_ye.dat" );

  eas_E_ye << "# T(MeV)="<<T*pc::k_MeV << " rho(g/cm^3)="<<rho << " ye_min="<<nulib_get_Yemin() << " ye_max="<<nulib_get_Yemax() << " Nye=" << ye_grid.size() << endl;
  eas_E_ye << setw(25) << "# ye" << setw(25) << "E(MeV)" << setw(25) << "emis(erg/cm^3/s/ster/Hz)" << setw(25) << "absopac(cm^2/g)" << setw(25) << "scatopac(cm^2/g)" << endl;

  for(int i=0; i<ye_grid.size(); i++){
    nulib_get_eas_arrays(rho, T, ye_grid[i], nulibID, emis, absopac, scatopac);
    nulib_get_pure_emis (rho_grid[i], T, ye, nulibID, pure_emis);
    for(int j=0; j<nu_grid.size(); j++){
      e = pure_emis[j];
      a = nu_grid.value_at(nu_grid[j], absopac);
      s = nu_grid.value_at(nu_grid[j], scatopac);
      eas_E_ye << setw(25) << ye_grid[i] <<setw(25) << nu_grid[j]*pc::h_MeV <<setw(25) << e  <<setw(25) << a/rho  <<setw(25) << s/rho << endl;
    }
    eas_E_ye << endl;
  }

  eas_E_ye.close();

  return 0;
}
