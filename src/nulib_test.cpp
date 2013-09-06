#include <iomanip>
#include <vector>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string.h>
#include <sstream>
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
  double rho      = atof(argv[2]);           // g/cm^3
  double T        = atof(argv[3])/pc::k_MeV; // K
  double ye       = atof(argv[4]);
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

  vector<double> ye_grid;
  nulib_get_Ye_array(ye_grid);
  grid_file << "Ye grid" << endl;
  for(int i=0; i<ye_grid.size(); i++) grid_file << ye_grid[i] << endl;
  grid_file << endl;
  
  vector<double> T_grid; // K
  nulib_get_T_array(T_grid);
  grid_file << "T grid (MeV)" << endl;
  for(int i=0; i<T_grid.size(); i++) grid_file << T_grid[i]*pc::k_MeV << endl;
  grid_file << endl;

  vector<double> rho_grid; // g/cm^3
  nulib_get_rho_array(rho_grid);
  grid_file << "rho grid (g/cm^3)" << endl;
  for(int i=0; i<rho_grid.size(); i++) grid_file << rho_grid[i] << endl;
  grid_file << endl;

  grid_file.close();

  // read in the number of species and groups in the table
  int ns = nulib_get_nspecies();
  int ng = nu_grid.size();
  cout << ns << " species and " << ng << " groups" << endl;

  //make vectors of appropriate sizes
  cout << "making vectors" << endl;
  vector<double> absopac (ng,0); // cm^-1
  vector<double> scatopac(ng,0); // cm^-1
  cdf_array emis;                // erg/cm^3/s
  emis.resize(ng);

  //===================//
  // SINGLE LINE PLOTS //
  //===================//
  double e,a,s;

  ofstream eas_rho;
  eas_rho.open("eas_rho.dat");
  cout << "generating eas_rho.dat" << endl;
  eas_rho << "# T(MeV):" << T*pc::k_MeV << " ye:" << ye << " E(MeV):" << myenergy << endl; 
  eas_rho << setw(20) << "# rho(g/cm^3)" << setw(20) << "emis(integrated)(erg/cm^3/s)" << setw(20) << "absopac(cm^-1)"<< setw(20) <<"scatopac(cm^-1)" << endl;
  for(int j=0; j<ng; j++){
    nulib_get_eas_arrays(rho_grid[j], T, ye, nulibID, emis, absopac, scatopac);
    e = emis.get(nu_grid.size()-1);
    a = nu_grid.value_at(myfreq, absopac);
    s = nu_grid.value_at(myfreq, scatopac);
    eas_rho << setw(20) << rho_grid[j] << setw(20) << e << setw(20) << a << setw(20) << s << endl;
  }

  ofstream eas_T;
  eas_T.open("eas_T.dat");
  cout << "generating eas_T.dat" << endl;
  eas_T << "# rho(g/cm^3):" << rho << " ye:" << ye << " E(MeV):" << myenergy << endl; 
  eas_T << setw(20) << "# T(MeV)" << setw(20) << "emis(integrated)(erg/cm^3/s)" << setw(20) << "absopac(cm^-1)" <<setw(20) << "scatopac(cm^-1)" << endl;
  for(int j=0; j<ng; j++){
    nulib_get_eas_arrays(rho, T_grid[j], ye, nulibID, emis, absopac, scatopac);
    e = emis.get(nu_grid.size()-1);
    a = nu_grid.value_at(myfreq, absopac);
    s = nu_grid.value_at(myfreq, scatopac);
    eas_T <<setw(20) << T_grid[j]*pc::k_MeV <<setw(20) << e <<setw(20) << a <<setw(20) << s << endl;
  }

  ofstream eas_ye;
  eas_ye.open("eas_ye.dat");
  cout << "generating eas_ye.dat" << endl;
  eas_ye << "# rho(g/cm^3):" << rho << " T(MeV):" << T*pc::k_MeV << " E(MeV):" << myenergy << endl; 
  eas_ye << setw(20) << "# ye" << setw(20) << "emis(integrated)(erg/cm^3/s)" << setw(20) << "absopac(cm^-1)" << setw(20) << "scatopac(cm^-1)" << endl;
  for(int j=0; j<ng; j++){
    nulib_get_eas_arrays(rho, T, ye_grid[j], nulibID, emis, absopac, scatopac);
    e = emis.get(nu_grid.size()-1);
    a = nu_grid.value_at(myfreq, absopac);
    s = nu_grid.value_at(myfreq, scatopac);
    eas_ye <<setw(20) << ye_grid[j] <<setw(20) << e <<setw(20) << a <<setw(20) << s << endl;
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

  eas_E_rho << "# T(MeV)="<<T*pc::k_MeV << " ye="<<ye << " rho_min(g/cm^3)="<<nulib_get_rhomin() << " rho_max(g/cm^3)="<<nulib_get_rhomax() << endl;
  eas_E_rho << setw(20) << "# rho(g/cm^3)" << setw(20) << "E(MeV)" << setw(20) << "emis(erg/cm^3/s/ster/Hz)" << setw(20) << "absopac(cm^-1)" << setw(20) << "scatopac(cm^-1)" << endl;

  for(int i=0; i<rho_grid.size(); i++){
    nulib_get_eas_arrays(rho_grid[i], T, ye, nulibID, emis, absopac, scatopac);
    for(int j=0; j<nu_grid.size(); j++){
      e = emis.get_value(j) / nu_grid.delta(j) / 4.0*pc::pi;
      a = nu_grid.value_at(nu_grid[j], absopac);
      s = nu_grid.value_at(nu_grid[j], scatopac);
      eas_E_rho << setw(20) << rho_grid[i] <<setw(20) << nu_grid[j]*pc::h_MeV <<setw(20) << e  <<setw(20) << a  <<setw(20) << s << endl;
    }
    eas_E_rho << endl;
  }

  eas_E_rho.close();


  // variation with T
  cout << "generating eas_E_T.dat" << endl;
  ofstream eas_E_T ("eas_E_T.dat" );

  eas_E_T << "# rho(g/cm^3)="<<rho << " ye="<<ye << " T_min(MeV)="<<nulib_get_Tmin()*pc::k_MeV << " T_max(MeV)="<<nulib_get_Tmax()*pc::k_MeV << endl;
  eas_E_T << setw(20) << "# T(MeV)" << setw(20) << "E(MeV)" << setw(20) << "emis(erg/cm^3/s/ster/Hz)" << setw(20) << "absopac(cm^-1)" << setw(20) << "scatopac(cm^-1)" << endl;

  for(int i=0; i<T_grid.size(); i++){
    nulib_get_eas_arrays(rho, T, ye, nulibID, emis, absopac, scatopac);
    for(int j=0; j<nu_grid.size(); j++){
      e = emis.get_value(j) / nu_grid.delta(j) / 4.0*pc::pi;
      a = nu_grid.value_at(nu_grid[j], absopac);
      s = nu_grid.value_at(nu_grid[j], scatopac);
      eas_E_T << setw(20) << T_grid[i]*pc::k_MeV << setw(20) << nu_grid[j]*pc::h_MeV << setw(20) << e << setw(20) << a << setw(20) << s << endl;
    }
    eas_E_T << endl;
  }

  eas_E_T.close();


  // variation with ye
  cout << "generating eas_E_ye.dat" << endl;
  ofstream eas_E_ye ("eas_E_ye.dat" );

  eas_E_ye << "# T(MeV)="<<T*pc::k_MeV << " rho(g/cm^3)="<<rho << " ye_min="<<nulib_get_Yemin << " ye_max="<<nulib_get_Yemax() << endl;
  eas_E_ye << setw(20) << "# ye" << setw(20) << "E(MeV)" << setw(20) << "emis(erg/cm^3/s/ster/Hz)" << setw(20) << "absopac(cm^-1)" << setw(20) << "scatopac(cm^-1)" << endl;

  for(int i=0; i<ye_grid.size(); i++){
    nulib_get_eas_arrays(rho, T, ye_grid[i], nulibID, emis, absopac, scatopac);
    for(int j=0; j<nu_grid.size(); j++){
      e = emis.get_value(j) / nu_grid.delta(j) / 4.0*pc::pi;
      a = nu_grid.value_at(nu_grid[j], absopac);
      s = nu_grid.value_at(nu_grid[j], scatopac);
      eas_E_ye << setw(20) << ye_grid[i] <<setw(20) << nu_grid[j]*pc::h_MeV <<setw(20) << e  <<setw(20) << a  <<setw(20) << s << endl;
    }
    eas_E_ye << endl;
  }

  eas_E_ye.close();


  ///////////////////

  // ofstream absopac_rho_ye("absopac_rho_ye.dat");
  // ofstream absopac_rho_T("absopac_rho_T.dat");
  // ofstream scatopac_rho_ye("scatopac_rho_ye.dat");
  // ofstream scatopac_rho_T("scatopac_rho_T.dat");
  // ofstream emis_T_ye("emis_T_ye.dat");
  // ofstream emis_T_rho("emis_T_rho.dat");

  // absopac_rho_ye  << "# [cm^-1] T(MeV):"           << T << endl;
  // absopac_rho_T   << "# [cm^-1] ye:"               << ye   << endl;
  // scatopac_rho_ye << "# [cm^-1] T(MeV):"           << T << endl;
  // scatopac_rho_T  << "# [cm^-1] ye:"               << ye   << endl;
  // emis_T_ye       << "# [erg/cm^3/s] rho(g/cm^3):" << rho  << endl; 
  // emis_T_rho      << "# [erg/cm^3/s] ye:"          << ye   << endl; 

  // double nu_min = 0;
  // double nu_max = nu_grid.size()-1;
  // double nu_med = nu_grid.size()/2;

  // absopac_rho_T   << "T(MeV)      rho(g/cm^3)  E(MeV)="<<nu_grid[nu_min]*pc::h_MeV<<" E(MeV)="<<nu_grid[nu_med]*pc::h_MeV<<" E(MeV)="<<nu_grid[nu_max]*pc::h_MeV << endl;
  // scatopac_rho_T  << "T(MeV)      rho(g/cm^3)  E(MeV)="<<nu_grid[nu_min]*pc::h_MeV<<" E(MeV)="<<nu_grid[nu_med]*pc::h_MeV<<" E(MeV)="<<nu_grid[nu_max]*pc::h_MeV << endl;
  // absopac_rho_ye  << "ye          rho(g/cm^3)  E(MeV)="<<nu_grid[nu_min]*pc::h_MeV<<" E(MeV)="<<nu_grid[nu_med]*pc::h_MeV<<" E(MeV)="<<nu_grid[nu_max]*pc::h_MeV << endl; 
  // scatopac_rho_ye << "ye          rho(g/cm^3)  E(MeV)="<<nu_grid[nu_min]*pc::h_MeV<<" E(MeV)="<<nu_grid[nu_med]*pc::h_MeV<<" E(MeV)="<<nu_grid[nu_max]*pc::h_MeV << endl;
  // emis_T_ye       << "ye          T(MeV)       (integrated over energies)" << endl; 
  // emis_T_rho      << "rho(g/cm^3) T(MeV)       (integrated over energies)" << endl;

  // // fill *_rho_T
  // for(int i=0; i<T_grid.size(); i++){
  //   for(int j=0; j<rho_grid.size(); j++){
  //     scatopac_rho_T << T_grid[i]*pc::k_MeV << " " << rho_grid[j] << " " << scatopac[nu_min] << " " << scatopac[nu_med] << " " << scatopac[nu_max] << endl;
  //   }
  //   absopac_rho_T   << endl;
  //   scatopac_rho_T  << endl;
  // }

  // // fill *_rho_ye
  // for(int i=0; i<ye_grid.size(); i++){
  //   for(int j=0; j<rho_grid.size(); j++){
  //     nulib_get_eas_arrays(rho_grid[j], T_grid[i], ye, nulibID, emis, absopac, scatopac);
  //     absopac_rho_ye  << ye_grid[i] << " " << rho_grid[j] << " " << absopac[nu_min]  << " " << absopac[nu_med]  << " " << absopac[nu_max]  << endl;
  //     scatopac_rho_ye << ye_grid[i] << " " << rho_grid[j] << " " << scatopac[nu_min] << " " << scatopac[nu_med] << " " << scatopac[nu_max] << endl;
  //   }
  //   absopac_rho_ye  << endl;
  //   scatopac_rho_ye << endl;
  // }


  // // fill emis_T_ye
  // for(int i=0; i<ye_grid.size(); i++){
  //   for(int j=0; j<T_grid.size(); j++){
  //     emis_T_ye << ye_grid[i] << " " << T_grid[j]*pc::k_MeV << " " << emis.N << endl;
  //   }
  //   emis_T_ye << endl;
  // }

  // // fill emis_T_rho
  // for(int i=0; i<rho_grid.size(); i++){
  //   for(int j=0; j<T_grid.size(); j++){
  //     emis_T_rho << rho_grid[i] << " " << T_grid[j]pc::k_MeV << " " << emis.N << endl;
  //   }
  //   emis_T_rho << endl;
  // }

  // // close the files
  // absopac_rho_ye.close();
  // absopac_rho_T.close();
  // scatopac_rho_ye.close();
  // scatopac_rho_T.close();
  // emis_T_ye.close();
  // emis_T_rho.close();


  return 0;
}
