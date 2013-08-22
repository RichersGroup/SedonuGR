#include <vector>
#include <iostream>
#include <stdlib.h>
#include <string.h>
#include "zone.h"
#include "nulib_interface.h"
#include "physical_constants.h"

namespace pc = physical_constants;

int main(){
  using namespace std;
  //set test inputs
  real rho = 1.0e10;
  real temp = 10.0/(1e-6*pc::k_ev);
  real ye = 0.35;
  real myenergy = 1.0;
  
  //read in the nulib table
  cout << "initializing nulib" << endl;
  string filename = "../external/NuLib/NuLib_LS220.h5";
  nulib_init(filename);
  
  // get the frequency array
  cout << "printing the frequency array" << endl;
  locate_array nu_grid;
  nulib_get_nu_grid(nu_grid);
  nu_grid.print();

  int ns = nulib_get_nspecies();
  int ng = nu_grid.size();
  cout << ns << " species and " << ng << "groups" << endl;

  //make vectors of appropriate sizes
  cout << "making vectors" << endl;
  vector<real> nut_absopac  (ng,0);
  vector<real> nut_scattopac(ng,0);
  cdf_array nut_emiss;
  nut_emiss.resize(ng);
  
  //fill the arrays
  cout << "filling array" << endl;
  int nulibID = 0;
  nulib_get_eas_arrays(rho, temp, ye, nulibID, nut_emiss, nut_absopac, nut_scattopac);
  
  //print the results
  cout << "rho: " << rho << endl;
  cout << "temp: " << temp << endl;
  cout << "ye: " << ye << endl;
  cout << "species group emissivity absopacity scatopacity" << endl;

  for(int j=0; j<ng; j++){
    cout << j << " " << nu_grid.x[j] << " " << " " << nut_emiss.get(j) << " " << nut_absopac[j] << " " << nut_scattopac[j] << endl;
  }

  return 0;
}
