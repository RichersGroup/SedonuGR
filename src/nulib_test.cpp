#include <vector>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string.h>
#include "zone.h"
#include "nulib_interface.h"
#include "physical_constants.h"

namespace pc = physical_constants;

int main(int argc, char* argv[]){
  using namespace std;
  //set test inputs
  double rho = 1.0e10;
  double temp = 10.0/(1e-6*pc::k_ev);
  double ye = 0.35;
  double myenergy = 1.0;
  int nulibID = 0;
  
  //read in the nulib table
  if(argc!=2){
    cout << "Usage: nulib_test.exe path_to_nulib_table.h5" << endl;
    exit(1);
  }
  cout << "initializing nulib" << endl;
  string filename = argv[1];
  nulib_init(filename);
  
  // get the frequency array
  locate_array nu_grid;
  nulib_get_nu_grid(nu_grid);

  // get the Ye array
  vector<double> ye_grid;
  nulib_get_Ye_array(ye_grid);
  
  // get the Temp array
  vector<double> T_grid;
  nulib_get_T_array(T_grid);

  // get the rho array
  vector<double> rho_grid;
  nulib_get_rho_array(rho_grid);

  // read in the number of species and groups in the table
  int ns = nulib_get_nspecies();
  int ng = nu_grid.size();
  cout << ns << " species and " << ng << " groups" << endl;

  //make vectors of appropriate sizes
  cout << "making vectors" << endl;
  vector<double> nut_absopac  (ng,0);
  vector<double> nut_scatopac(ng,0);
  cdf_array nut_emis;
  nut_emis.resize(ng);

  // nu dependence
  ofstream eas_nu;
  eas_nu.open("eas_nu.dat");
  cout << "Now generating *_nu.dat" << endl;
  eas_nu << "# rho:" << rho << " T:" << temp << " ye:" << ye << endl; 
  eas_nu << "# index - nu - emis - absopac - scatopac" << endl;
  nulib_get_eas_arrays(rho, temp, ye, nulibID, nut_emis, nut_absopac, nut_scatopac);
  for(int j=0; j<ng; j++){
    eas_nu << j << " " << nu_grid.x[j] << " " << nut_emis.get(j) << " " << nut_absopac[j] << " " << nut_scatopac[j] << endl;
  }

  // prepare other files
  cout << "Now generating other files." << endl;
  ofstream absopac_rho_ye("absopac_rho_ye.dat");
  ofstream absopac_rho_T("absopac_rho_T.dat");
  ofstream scatopac_rho_ye("scatopac_rho_ye.dat");
  ofstream scatopac_rho_T("scatopac_rho_T.dat");
  ofstream emis_rho_ye("emis_rho_ye.dat");
  ofstream emis_rho_T("emis_rho_T.dat");

  absopac_rho_ye  << "# rho:" << rho << " T:" << temp << " ye:" << ye << endl; 
  absopac_rho_T   << "# rho:" << rho << " T:" << temp << " ye:" << ye << endl; 
  scatopac_rho_ye << "# rho:" << rho << " T:" << temp << " ye:" << ye << endl; 
  scatopac_rho_T  << "# rho:" << rho << " T:" << temp << " ye:" << ye << endl; 
  emis_rho_ye     << "# rho:" << rho << " T:" << temp << " ye:" << ye << endl; 
  emis_rho_T      << "# rho:" << rho << " T:" << temp << " ye:" << ye << endl; 

  double nu_min = 0;
  double nu_max = nu_grid.size()-1;
  double nu_med = nu_grid.size()/2;
  double T_med  =  T_grid.size()/2;
  double ye_med = ye_grid.size()/2;

  absopac_rho_ye  << "rho  ye  nu="<<nu_grid[nu_min]<<" nu="<<nu_grid[nu_med]<<" nu="<<nu_grid[nu_max] << endl; 
  absopac_rho_T   << "rho  T   nu="<<nu_grid[nu_min]<<" nu="<<nu_grid[nu_med]<<" nu="<<nu_grid[nu_max] << endl;
  scatopac_rho_ye << "rho  ye  nu="<<nu_grid[nu_min]<<" nu="<<nu_grid[nu_med]<<" nu="<<nu_grid[nu_max] << endl;
  scatopac_rho_T  << "rho  T   nu="<<nu_grid[nu_min]<<" nu="<<nu_grid[nu_med]<<" nu="<<nu_grid[nu_max] << endl;
  emis_rho_ye     << "rho  ye  nu="<<nu_grid[nu_min]<<" nu="<<nu_grid[nu_med]<<" nu="<<nu_grid[nu_max] << endl; 
  emis_rho_T      << "rho  T   nu="<<nu_grid[nu_min]<<" nu="<<nu_grid[nu_med]<<" nu="<<nu_grid[nu_max] << endl;

  // fill other files
  for(int i=0; i<rho_grid.size(); i++)
  {
    for(int j=0; j<ye_grid.size(); j++)
    {
      nulib_get_eas_arrays(rho_grid[i], T_grid[T_med], ye_grid[j], nulibID, nut_emis, nut_absopac, nut_scatopac);
      absopac_rho_ye  << rho_grid[i] << " " << ye_grid[j] << " " << nut_absopac[nu_min]  << " " << nut_absopac[nu_med]  << " " << nut_absopac[nu_max]  << endl;
      scatopac_rho_ye << rho_grid[i] << " " << ye_grid[j] << " " << nut_scatopac[nu_min] << " " << nut_scatopac[nu_med] << " " << nut_scatopac[nu_max] << endl;
      emis_rho_ye     << rho_grid[i] << " " << ye_grid[j] << " " << nut_emis.get(nu_min) << " " << nut_emis.get(nu_med) << " " << nut_emis.get(nu_max) << endl;
    }
    for(int j=0; j<T_grid.size(); j++)
    {
      nulib_get_eas_arrays(rho_grid[i], T_grid[j], ye_grid[ye_med], nulibID, nut_emis, nut_absopac, nut_scatopac);
      absopac_rho_T  << rho_grid[i] << " " << T_grid[j] << " " << nut_absopac[nu_min]  << " " << nut_absopac[nu_med]  << " " << nut_absopac[nu_max]  << endl;
      scatopac_rho_T << rho_grid[i] << " " << T_grid[j] << " " << nut_scatopac[nu_min] << " " << nut_scatopac[nu_med] << " " << nut_scatopac[nu_max] << endl;
      emis_rho_T     << rho_grid[i] << " " << T_grid[j] << " " << nut_emis.get(nu_min) << " " << nut_emis.get(nu_med) << " " << nut_emis.get(nu_max) << endl;
    }

    // space makes gnuplot interpret correctly
    absopac_rho_ye  << endl;
    scatopac_rho_ye << endl;
    emis_rho_ye     << endl;
    absopac_rho_T   << endl;
    scatopac_rho_T  << endl;
    emis_rho_T      << endl;
  }

  // close the files
  absopac_rho_ye.close();
  absopac_rho_T.close();
  scatopac_rho_ye.close();
  scatopac_rho_T.close();
  emis_rho_ye.close();
  emis_rho_T.close();
  eas_nu.close();

  return 0;
}
