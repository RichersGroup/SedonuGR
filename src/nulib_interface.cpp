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
// The global variables that will be used, independent of fortran compiler version
int     nulibtable_number_species;
int     nulibtable_number_easvariables;
int     nulibtable_number_groups;
int     nulibtable_nrho;
int     nulibtable_ntemp;
int     nulibtable_nye;
double* nulibtable_energies;
double* nulibtable_ewidths;
double* nulibtable_ebottom;
double* nulibtable_logrho;
double* nulibtable_logtemp;
double* nulibtable_ye;
double  nulibtable_logtemp_min;
double  nulibtable_logtemp_max;
double  nulibtable_logrho_min;
double  nulibtable_logrho_max;
double  nulibtable_ye_min;
double  nulibtable_ye_max;

// The format of the fortran variables the fortran compiler provides
// assumes C and Fortran compilers are the same
// To be copied into the universal globals if intel compiler
#ifdef __INTEL_COMPILER
extern int     nulibtable_mp_nulibtable_number_species_;
extern int     nulibtable_mp_nulibtable_number_easvariables_;
extern int     nulibtable_mp_nulibtable_number_groups_;
extern int     nulibtable_mp_nulibtable_nrho_;
extern int     nulibtable_mp_nulibtable_ntemp_;
extern int     nulibtable_mp_nulibtable_nye_;
extern double* nulibtable_mp_nulibtable_energies_;
extern double* nulibtable_mp_nulibtable_ewidths_;
extern double* nulibtable_mp_nulibtable_ebottom_;
extern double* nulibtable_mp_nulibtable_logrho_;
extern double* nulibtable_mp_nulibtable_logtemp_;
extern double* nulibtable_mp_nulibtable_ye_;
extern double  nulibtable_mp_nulibtable_logtemp_min_;
extern double  nulibtable_mp_nulibtable_logtemp_max_;
extern double  nulibtable_mp_nulibtable_logrho_min_;
extern double  nulibtable_mp_nulibtable_logrho_max_;
extern double  nulibtable_mp_nulibtable_ye_min_;
extern double  nulibtable_mp_nulibtable_ye_max_;
#elif defined __GNUC__
extern int     __nulibtable_MOD_nulibtable_number_species;
extern int     __nulibtable_MOD_nulibtable_number_easvariables;
extern int     __nulibtable_MOD_nulibtable_number_groups;
extern int     __nulibtable_MOD_nulibtable_nrho;
extern int     __nulibtable_MOD_nulibtable_ntemp;
extern int     __nulibtable_MOD_nulibtable_nye;
extern double* __nulibtable_MOD_nulibtable_energies;
extern double* __nulibtable_MOD_nulibtable_ewidths;
extern double* __nulibtable_MOD_nulibtable_ebottom;
extern double* __nulibtable_MOD_nulibtable_logrho;
extern double* __nulibtable_MOD_nulibtable_logtemp;
extern double* __nulibtable_MOD_nulibtable_ye;
extern double  __nulibtable_MOD_nulibtable_logtemp_min;
extern double  __nulibtable_MOD_nulibtable_logtemp_max;
extern double  __nulibtable_MOD_nulibtable_logrho_min;
extern double  __nulibtable_MOD_nulibtable_logrho_max;
extern double  __nulibtable_MOD_nulibtable_ye_min;
extern double  __nulibtable_MOD_nulibtable_ye_max;
#elif defined __PGI
#else
#warning "The fortran interface is only configured for Intel and GNU compilers. Attempting default variable names. If this does not work you must modify src/nulib_interface.cpp to get your C++ and Fortran compilers to play nicely together."
#endif


/**************************************/
/* set the universal global variables */
/**************************************/
void set_globals(){
#ifdef __INTEL_COMPILER
  nulibtable_number_species       = nulibtable_mp_nulibtable_number_species_;
  nulibtable_number_easvariables  = nulibtable_mp_nulibtable_number_easvariables_;
  nulibtable_number_groups        = nulibtable_mp_nulibtable_number_groups_;
  nulibtable_nrho                 = nulibtable_mp_nulibtable_nrho_;
  nulibtable_ntemp                = nulibtable_mp_nulibtable_ntemp_;
  nulibtable_nye                  = nulibtable_mp_nulibtable_nye_;
  nulibtable_energies             = nulibtable_mp_nulibtable_energies_;
  nulibtable_ewidths              = nulibtable_mp_nulibtable_ewidths_;
  nulibtable_ebottom              = nulibtable_mp_nulibtable_ebottom_;
  nulibtable_logrho               = nulibtable_mp_nulibtable_logrho_;
  nulibtable_logtemp              = nulibtable_mp_nulibtable_logtemp_;
  nulibtable_ye                   = nulibtable_mp_nulibtable_ye_;
  nulibtable_logtemp_min          = nulibtable_mp_nulibtable_logtemp_min_;
  nulibtable_logtemp_max          = nulibtable_mp_nulibtable_logtemp_max_;
  nulibtable_logrho_min           = nulibtable_mp_nulibtable_logrho_min_;
  nulibtable_logrho_max           = nulibtable_mp_nulibtable_logrho_max_;
  nulibtable_ye_min               = nulibtable_mp_nulibtable_ye_min_;
  nulibtable_ye_max               = nulibtable_mp_nulibtable_ye_max_;
#elif defined __GNUC__
  nulibtable_number_species      = __nulibtable_MOD_nulibtable_number_species;
  nulibtable_number_easvariables = __nulibtable_MOD_nulibtable_number_easvariables;
  nulibtable_number_groups       = __nulibtable_MOD_nulibtable_number_groups;
  nulibtable_nrho                = __nulibtable_MOD_nulibtable_nrho;
  nulibtable_ntemp               = __nulibtable_MOD_nulibtable_ntemp;
  nulibtable_nye                 = __nulibtable_MOD_nulibtable_nye;
  nulibtable_energies            = __nulibtable_MOD_nulibtable_energies;
  nulibtable_ewidths             = __nulibtable_MOD_nulibtable_ewidths;
  nulibtable_ebottom             = __nulibtable_MOD_nulibtable_ebottom;
  nulibtable_logrho              = __nulibtable_MOD_nulibtable_logrho;
  nulibtable_logtemp             = __nulibtable_MOD_nulibtable_logtemp;
  nulibtable_ye                  = __nulibtable_MOD_nulibtable_ye;
  nulibtable_logtemp_min         = __nulibtable_MOD_nulibtable_logtemp_min;
  nulibtable_logtemp_max         = __nulibtable_MOD_nulibtable_logtemp_max;
  nulibtable_logrho_min          = __nulibtable_MOD_nulibtable_logrho_min;
  nulibtable_logrho_max          = __nulibtable_MOD_nulibtable_logrho_max;
  nulibtable_ye_min              = __nulibtable_MOD_nulibtable_ye_min;
  nulibtable_ye_max              = __nulibtable_MOD_nulibtable_ye_max;
#endif
}

/**********************/
/* nulib_get_nspecies */
/**********************/
int nulib_get_nspecies(){
  return nulibtable_number_species;
}

/**************/
/* nulib_init */
/**************/
void nulib_init(string filename){
  nulibtable_reader_((char*)filename.c_str(), filename.length());
  set_globals();
}


/**********************/
/* get_nut_eas_arrays */
/**********************/
// leave serial - called per grid zone
void nulib_get_eas_arrays(real rho,                     // g/cm^3 
			  real temp,                    // K
			  real ye, int nulibID,
			  cdf_array& nut_emiss,         // erg/cm^3/s/ster
			  vector<real>& nut_absopac,    // cm^-1   
			  vector<real>& nut_scatopac){  // cm^-1
  
  int nvars    = nulibtable_number_easvariables;
  int ngroups  = nulibtable_number_groups;
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
  if(log10(rho) < nulibtable_logrho_min) 
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
  widths.assign(nulibtable_ewidths, 
		nulibtable_ewidths + nulibtable_number_groups);
 
  int nvars    = nulibtable_number_easvariables;
  int ngroups  = nulibtable_number_groups;
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
  if(log10(rho) < nulibtable_logrho_min) 
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
  nu_grid.x.assign(nulibtable_energies, 
		   nulibtable_energies + nulibtable_number_groups);

  // convert from MeV to frequency using the Planck constant
  for(int i=0; i<nu_grid.size(); i++) nu_grid.x[i] /= pc::h_MeV;
  nu_grid.do_log_interpolate = 1;
  nu_grid.min = nulibtable_ebottom[0];
}


/*****************************/
/* get the rho, T, Ye arrays */
/*****************************/
void nulib_get_rho_array(vector<double>& array){ // g/cm^3
  array.assign(nulibtable_logrho,
	       nulibtable_logrho  + nulibtable_nrho);
  for(int i=0; i<array.size(); i++) array[i] = pow(10.0, array[i]);
}
void nulib_get_T_array(vector<double>& array){ // K
  array.assign(nulibtable_logtemp,
	       nulibtable_logtemp + nulibtable_ntemp);
  for(int i=0; i<array.size(); i++) array[i] = pow(10.0, array[i]) / pc::k_MeV; // convert from MeV to K
}
void nulib_get_Ye_array(vector<double>& array){
  array.assign(nulibtable_ye,
	       nulibtable_ye      + nulibtable_nye);
}

/*************************/
/* nulib_get_{*min,*max} */
/*************************/
double nulib_get_Tmin()   {return pow(10,nulibtable_logtemp_min) / pc::k_MeV;} //convert from MeV to K
double nulib_get_Tmax()   {return pow(10,nulibtable_logtemp_max) / pc::k_MeV;} //convert from MeV to K
double nulib_get_rhomin() {return pow(10,nulibtable_logrho_min);}
double nulib_get_rhomax() {return pow(10,nulibtable_logrho_max);}
double nulib_get_Yemin()  {return nulibtable_ye_min;}
double nulib_get_Yemax()  {return nulibtable_ye_max;}
