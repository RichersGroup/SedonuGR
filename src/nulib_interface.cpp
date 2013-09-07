#include <vector>
#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <math.h>
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
extern int __nulibtable_MOD_nulibtable_nrho;
extern int __nulibtable_MOD_nulibtable_ntemp;
extern int __nulibtable_MOD_nulibtable_nye;
extern double* __nulibtable_MOD_nulibtable_energies;
extern double* __nulibtable_MOD_nulibtable_ewidths;
extern double* __nulibtable_MOD_nulibtable_logrho;
extern double* __nulibtable_MOD_nulibtable_logtemp;
extern double* __nulibtable_MOD_nulibtable_ye;
extern double __nulibtable_MOD_nulibtable_logtemp_min;
extern double __nulibtable_MOD_nulibtable_logtemp_max;
extern double __nulibtable_MOD_nulibtable_logrho_min;
extern double __nulibtable_MOD_nulibtable_logrho_max;
extern double __nulibtable_MOD_nulibtable_ye_min;
extern double __nulibtable_MOD_nulibtable_ye_max;


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
void nulib_get_eas_arrays(real rho,                     // g/cm^3 
			  real temp,                    // K
			  real ye, int nulibID,
			  cdf_array& nut_emiss,         // erg/cm^3/s/ster
			  vector<real>& nut_absopac,    // cm^-1   
			  vector<real>& nut_scatopac){  // cm^-1
  
  int nvars    = __nulibtable_MOD_nulibtable_number_easvariables;
  int ngroups  = __nulibtable_MOD_nulibtable_number_groups;
  // apparently it's valid to declare array sizes at runtime like this... if it breaks, use malloc
  double eas_energy[nvars][ngroups]; //[0][*] = energy-bin-integrated emissivity (erg/cm^3/s/ster). 
                                     //[1][*] = absorption opacity (1/cm)
                                     //[2][*] = scattering opacity (1/cm)

  //check sizes match and we are within the table boundaries
  if(nvars != 3){
    cout << "ERROR: nulibtable has different number of eas variables than 3" << endl;
    exit(EXIT_FAILURE);
  }
  if(ngroups != nut_absopac.size()){
    cout << "ERROR: nulibtable has different number of energies than the eas arrays." << endl;
    exit(EXIT_FAILURE);
  }

  // fetch the relevant table from nulib. NuLib only accepts doubles.
  double temp_MeV = temp * pc::k_MeV; // MeV
  int lns = nulibID+1;                // fortran array indices start with 1
  
  // If the density is too low, just set everything to zero
  if(log10(rho) < __nulibtable_MOD_nulibtable_logrho_min) 
    for(int j=0; j<ngroups; j++){      
      nut_emiss.set_value(j, 0);
      nut_absopac [j] =      0;
      nut_scatopac[j] =      0;
    }
  
  // Otherwise, fill with the appropriate values
  else{
    // nulib's emis is bin integrated. Summing these values would correspond to our intended nut_emiss value at the bin top.
    // must rebin to get the integrated value to be at the same location as the opacities. (CDF value corresponds to emission rate at or below that energy)
    nulibtable_single_species_range_energy_(&rho, &temp_MeV, &ye, &lns,
					    (double*)eas_energy, &ngroups, &nvars);
    double last_emis = 0;
    double this_emis = 0;
    for(int j=0; j<ngroups; j++){
      last_emis = this_emis;
      this_emis = eas_energy[0][j];
      nut_emiss.set_value(j, 0.5*(last_emis + this_emis) );
      nut_absopac [j] =      eas_energy[1][j];
      nut_scatopac[j] =      eas_energy[2][j];
    }
  }
  nut_emiss.N = 1.0;
}


/***********************/
/* nulib_get_pure_emis */
/***********************/
void nulib_get_pure_emis(real rho,                     // g/cm^3 
			 real temp,                    // K
			 real ye, int nulibID,
			 vector<double>& nut_emiss){   // erg/cm^3/s/ster/Hz
  
  vector<double> widths;
  widths.assign(__nulibtable_MOD_nulibtable_ewidths, 
		__nulibtable_MOD_nulibtable_ewidths + __nulibtable_MOD_nulibtable_number_groups);
 
  int nvars    = __nulibtable_MOD_nulibtable_number_easvariables;
  int ngroups  = __nulibtable_MOD_nulibtable_number_groups;
  // apparently it's valid to declare array sizes at runtime like this... if it breaks, use malloc
  double eas_energy[nvars][ngroups]; //[0][*] = energy-bin-integrated emissivity (erg/cm^3/s/ster). 
                                     //[1][*] = absorption opacity (1/cm)
                                     //[2][*] = scattering opacity (1/cm)

  //check sizes match and we are within the table boundaries
  if(nvars != 3){
    cout << "ERROR: nulibtable has different number of eas variables than 3" << endl;
    exit(EXIT_FAILURE);
  }
  if(ngroups != nut_emiss.size()){
    cout << "ERROR: nulibtable has different number of energies than the eas arrays." << endl;
    exit(EXIT_FAILURE);
  }

  // fetch the relevant table from nulib. NuLib only accepts doubles.
  double temp_MeV = temp * pc::k_MeV; // MeV
  int lns = nulibID+1;                // fortran array indices start with 1
  
  // If the density is too low, just set everything to zero
  if(log10(rho) < __nulibtable_MOD_nulibtable_logrho_min) 
    for(int j=0; j<ngroups; j++) nut_emiss[j] = 0;
  
  // Otherwise, fill with the appropriate values
  else{
    nulibtable_single_species_range_energy_(&rho, &temp_MeV, &ye, &lns,
					    (double*)eas_energy, &ngroups, &nvars);
    for(int j=0; j<ngroups; j++) nut_emiss[j] = eas_energy[0][j] / (widths[j]/pc::h_MeV);
  }
}


/*********************/
/* nulib_get_nu_grid */
/*********************/
// Fill in the locate array with values of the array stored in the fortran module
void nulib_get_nu_grid(locate_array& nu_grid){ // Hz
  // assign values from the NuLib module to nu_grid
  nu_grid.x.assign(__nulibtable_MOD_nulibtable_energies, 
		   __nulibtable_MOD_nulibtable_energies + __nulibtable_MOD_nulibtable_number_groups);

  // convert from MeV to frequency using the Planck constant
  for(int i=0; i<nu_grid.size(); i++) nu_grid.x[i] /= pc::h_MeV;
  nu_grid.do_log_interpolate = 1;
}


/*****************************/
/* get the rho, T, Ye arrays */
/*****************************/
void nulib_get_rho_array(vector<double>& array){ // g/cm^3
  array.assign(__nulibtable_MOD_nulibtable_logrho,
	       __nulibtable_MOD_nulibtable_logrho  + __nulibtable_MOD_nulibtable_nrho);
  for(int i=0; i<array.size(); i++) array[i] = pow(10.0, array[i]);
}
void nulib_get_T_array(vector<double>& array){ // K
  array.assign(__nulibtable_MOD_nulibtable_logtemp,
	       __nulibtable_MOD_nulibtable_logtemp + __nulibtable_MOD_nulibtable_ntemp);
  for(int i=0; i<array.size(); i++) array[i] = pow(10.0, array[i]) / pc::k_MeV; // convert from MeV to K
}
void nulib_get_Ye_array(vector<double>& array){
  array.assign(__nulibtable_MOD_nulibtable_ye,
	       __nulibtable_MOD_nulibtable_ye      + __nulibtable_MOD_nulibtable_nye);
}

/*************************/
/* nulib_get_{*min,*max} */
/*************************/
double nulib_get_Tmin()   {return pow(10,__nulibtable_MOD_nulibtable_logtemp_min) / pc::k_MeV;} //convert from MeV to K
double nulib_get_Tmax()   {return pow(10,__nulibtable_MOD_nulibtable_logtemp_max) / pc::k_MeV;} //convert from MeV to K
double nulib_get_rhomin() {return pow(10,__nulibtable_MOD_nulibtable_logrho_min);}
double nulib_get_rhomax() {return pow(10,__nulibtable_MOD_nulibtable_logrho_max);}
double nulib_get_Yemin()  {return __nulibtable_MOD_nulibtable_ye_min;}
double nulib_get_Yemax()  {return __nulibtable_MOD_nulibtable_ye_max;}
