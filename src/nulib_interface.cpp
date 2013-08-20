#include <vector>
#include <iostream>
#include <stdlib.h>
#include <string.h>
#include "zone.h"
#include "nulib_interface.h"
#include "physical_constants.h"

using namespace std;
namespace pc = physical_constants;

// These are fortran functions and module variables in nulib.a
extern "C"{
  void nulibtable_range_species_range_energy_(double*, //rho
					      double*, //temp
					      double*, //ye
					      double*, //eas_species_energy (3D array)
					      int*,    //number of species (3,5,6)
					      int*,    //number of groups
					      int*);   //number of easvariables (3)

  void nulibtable_single_species_range_energy_(double*, //rho
					       double*, //temp
					       double*, //Ye
					       int*,    //species number
					       double*, //eas_energy (2D array)
					       int*,    //number of groups
					       int*);   //number of easvariables(3)
					       
  void nulibtable_reader_(char*,int);
}
// TODO - I think these names are compiler-dependent (this is gfortran). Check with ifort.
extern int __nulibtable_MOD_nulibtable_number_species;
extern int __nulibtable_MOD_nulibtable_number_easvariables;
extern int __nulibtable_MOD_nulibtable_number_groups;
extern double* __nulibtable_MOD_nulibtable_energies;

/**********************/
/* nulib_get_nspecies */
/**********************/
int nulib_get_nspecies(){
  return __nulibtable_MOD_nulibtable_number_species;
}

/**************/
/* nulib_init */
/**************/
void nulib_init(string filename){
  nulibtable_reader_((char*)filename.c_str(), filename.length());
}


/**********************/
/* get_nut_eas_arrays */
/**********************/
void nulib_get_eas_arrays(real rho, real temp, real ye, int nulibID,
			  cdf_array& nut_emiss, vector<real>& nut_absopac, vector<real>& nut_scatopac){
  
  int nvars    = __nulibtable_MOD_nulibtable_number_easvariables;
  int ngroups  = __nulibtable_MOD_nulibtable_number_groups;
  //  double eas_species_energy[nvars][ngroups][nspecies];
  // apparently it's valid to declare array sizes at runtime like this... if it breaks, use malloc
  double eas_energy[nvars][ngroups];

  //check sizes match and we are withing the table boundaries
  if(nvars != 3){
    cout << "ERROR: nulibtable has different number of eas variables than 3" << endl;
    exit(EXIT_FAILURE);
  }
  if(ngroups != nut_absopac.size()){
    cout << "ERROR: nulibtable has different number of energies than the eas arrays." << endl;
    exit(EXIT_FAILURE);
  }

  // fetch the relevant table from nulib. NuLib only accepts doubles.
  double rhotmp = rho;
  double temptmp = temp * 1e-6*pc::k_ev; // convert temperature to MeV
  double yetmp = ye;
  int lns = nulibID+1;                   // fortran array indices start with 1
  nulibtable_single_species_range_energy_(&rhotmp, &temptmp, &yetmp, &lns,
					 (double*)eas_energy, &ngroups, &nvars);

  //fill the vectors with appropriate values
  for(int j=0; j<ngroups; j++){
    nut_emiss.set_value(j, eas_energy[0][j]);
    nut_absopac [j] =      eas_energy[1][j];
    nut_scatopac[j] =      eas_energy[2][j];
  }
}


/*********************/
/* nulib_get_nu_grid */
/*********************/
// Fill in the locate array with values of the array stored in the fortran module
void nulib_get_nu_grid(locate_array& nu_grid){
  // assign values from the NuLib module to nu_grid
  nu_grid.x.assign(__nulibtable_MOD_nulibtable_energies, 
		   __nulibtable_MOD_nulibtable_energies + __nulibtable_MOD_nulibtable_number_groups);

  // convert from MeV to frequency using the Planck constant
  for(int i=0; i<nu_grid.size(); i++) nu_grid.x[i] /= pc::h;
}
