/*
//  Copyright (c) 2015, California Institute of Technology and the Regents
//  of the University of California, based on research sponsored by the
//  United States Department of Energy. All rights reserved.
//
//  This file is part of Sedonu.
//
//  Sedonu is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  Neither the name of the California Institute of Technology (Caltech)
//  nor the University of California nor the names of its contributors 
//  may be used to endorse or promote products derived from this software
//  without specific prior written permission.
//
//  Sedonu is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with Sedonu.  If not, see <http://www.gnu.org/licenses/>.
//
*/

#include "global_options.h"
#include <vector>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include "zone.h"
#include "nulib_interface.h"

// These are fortran functions and module variables in nulib.a
extern "C"{
void nulibtable_range_species_range_energy_(
		double*, //rho
		double*, //temp
		double*, //ye
		double*, //eas_species_energy (3D array)
		int*,    //number of species (3,5,6)
		int*,    //number of groups
		int*);   //number of easvariables (3)

void nulibtable_single_species_range_energy_(
		double*, //rho
		double*, //temp
		double*, //Ye
		int*,    //species number
		double*, //eas_energy (2D array)
		int*,    //number of groups
		int*);   //number of easvariables(3)
void nulibtable_epannihil_single_species_range_energy_(
		double* temp,  // MeV
		double* eta,
		int* lns,
		double* eas,
		int* ngroups1,
		int* ngroups2,
		int* n_phis);

void nulibtable_reader_(char*,int*,int*,int);

void read_eos_table_(char* filename);

void set_eos_variables_(double* eos_variables);
}
// The global variables that will be used, independent of fortran compiler version
int     nulibtable_number_species;
int     nulibtable_number_easvariables;
int     nulibtable_number_groups;
int     nulibtable_nrho;
int     nulibtable_ntemp;
int     nulibtable_nye;
int 	nulibtable_nItemp;
int 	nulibtable_nIeta;
double* nulibtable_energies;
double* nulibtable_ewidths;
double* nulibtable_ebottom;
double* nulibtable_etop;
double* nulibtable_logrho;
double* nulibtable_logtemp;
double* nulibtable_ye;
double* nulibtable_logItemp;
double* nulibtable_logIeta;
double  nulibtable_logtemp_min;
double  nulibtable_logtemp_max;
double  nulibtable_logrho_min;
double  nulibtable_logrho_max;
double  nulibtable_ye_min;
double  nulibtable_ye_max;
double  nulibtable_logItemp_min;
double  nulibtable_logItemp_max;
double  nulibtable_logIeta_min;
double  nulibtable_logIeta_max;
int     nulib_total_eos_variables;

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
extern int	   nulibtable_mp_nulibtable_nitemp_;
extern int     nulibtable_mp_nulibtable_nieta_;
extern double* nulibtable_mp_nulibtable_energies_;
extern double* nulibtable_mp_nulibtable_ewidths_;
extern double* nulibtable_mp_nulibtable_ebottom_;
extern double* nulibtable_mp_nulibtable_etop_;
extern double* nulibtable_mp_nulibtable_logrho_;
extern double* nulibtable_mp_nulibtable_logtemp_;
extern double* nulibtable_mp_nulibtable_ye_;
extern double* nulibtable_mp_nulibtable_logitemp;
extern double* nulibtable_mp_nulibtable_logieta;
extern double  nulibtable_mp_nulibtable_logtemp_min_;
extern double  nulibtable_mp_nulibtable_logtemp_max_;
extern double  nulibtable_mp_nulibtable_logrho_min_;
extern double  nulibtable_mp_nulibtable_logrho_max_;
extern double  nulibtable_mp_nulibtable_ye_min_;
extern double  nulibtable_mp_nulibtable_ye_max_;
extern double  nulibtable_mp_nulibtable_logitemp_min_;
extern double  nulibtable_mp_nulibtable_logitemp_max_;
extern double  nulibtable_mp_nulibtable_logieta_min_;
extern double  nulibtable_mp_nulibtable_logieta_max_;
extern int     nulib_mp_total_eos_variables_;
#elif defined __GNUC__
extern int     __nulibtable_MOD_nulibtable_number_species;
extern int     __nulibtable_MOD_nulibtable_number_easvariables;
extern int     __nulibtable_MOD_nulibtable_number_groups;
extern int     __nulibtable_MOD_nulibtable_nrho;
extern int     __nulibtable_MOD_nulibtable_ntemp;
extern int     __nulibtable_MOD_nulibtable_nye;
extern int	   __nulibtable_MOD_nulibtable_nitemp;
extern int     __nulibtable_MOD_nulibtable_nieta;
extern double* __nulibtable_MOD_nulibtable_energies;
extern double* __nulibtable_MOD_nulibtable_ewidths;
extern double* __nulibtable_MOD_nulibtable_ebottom;
extern double* __nulibtable_MOD_nulibtable_etop;
extern double* __nulibtable_MOD_nulibtable_logrho;
extern double* __nulibtable_MOD_nulibtable_logtemp;
extern double* __nulibtable_MOD_nulibtable_ye;
extern double* __nulibtable_MOD_nulibtable_logitemp;
extern double* __nulibtable_MOD_nulibtable_logieta;
extern double  __nulibtable_MOD_nulibtable_logtemp_min;
extern double  __nulibtable_MOD_nulibtable_logtemp_max;
extern double  __nulibtable_MOD_nulibtable_logrho_min;
extern double  __nulibtable_MOD_nulibtable_logrho_max;
extern double  __nulibtable_MOD_nulibtable_ye_min;
extern double  __nulibtable_MOD_nulibtable_ye_max;
extern double  __nulibtable_MOD_nulibtable_logitemp_min;
extern double  __nulibtable_MOD_nulibtable_logitemp_max;
extern double  __nulibtable_MOD_nulibtable_logieta_min;
extern double  __nulibtable_MOD_nulibtable_logieta_max;
extern int     __nulib_MOD_total_eos_variables;
#elif defined __PGI
#else
#warning "The fortran interface is only configured for Intel and GNU compilers. Attempting default variable names. If this does not work you must modify src/nulib_interface.cpp to get your C++ and Fortran compilers to play nicely together."
#endif


/**************************************/
/* set the universal global variables */
/**************************************/
void nulibtable_set_globals(){
#ifdef __INTEL_COMPILER
	nulibtable_number_species       = nulibtable_mp_nulibtable_number_species_;
	nulibtable_number_easvariables  = nulibtable_mp_nulibtable_number_easvariables_;
	nulibtable_number_groups        = nulibtable_mp_nulibtable_number_groups_;
	nulibtable_nrho                 = nulibtable_mp_nulibtable_nrho_;
	nulibtable_ntemp                = nulibtable_mp_nulibtable_ntemp_;
	nulibtable_nye                  = nulibtable_mp_nulibtable_nye_;
	nulibtable_nItemp				= nulibtable_mp_nulibtable_nitemp_;
	nulibtable_nIeta				= nulibtable_mp_nulibtable_nieta_;
	nulibtable_energies             = nulibtable_mp_nulibtable_energies_;
	nulibtable_ewidths              = nulibtable_mp_nulibtable_ewidths_;
	nulibtable_ebottom              = nulibtable_mp_nulibtable_ebottom_;
	nulibtable_etop                 = nulibtable_mp_nulibtable_etop_;
	nulibtable_logrho               = nulibtable_mp_nulibtable_logrho_;
	nulibtable_logtemp              = nulibtable_mp_nulibtable_logtemp_;
	nulibtable_ye                   = nulibtable_mp_nulibtable_ye_;
	nulibtable_logItemp             = nulibtable_mp_nulibtable_logitemp_;
	nulibtable_logIeta              = nulibtable_mp_nulibtable_logieta_;
	nulibtable_logtemp_min          = nulibtable_mp_nulibtable_logtemp_min_;
	nulibtable_logtemp_max          = nulibtable_mp_nulibtable_logtemp_max_;
	nulibtable_logrho_min           = nulibtable_mp_nulibtable_logrho_min_;
	nulibtable_logrho_max           = nulibtable_mp_nulibtable_logrho_max_;
	nulibtable_ye_min               = nulibtable_mp_nulibtable_ye_min_;
	nulibtable_ye_max               = nulibtable_mp_nulibtable_ye_max_;
	nulibtable_logItemp_min		    = nulibtable_mp_nulibtable_logitemp_min_;
	nulibtable_logItemp_max         = nulibtable_mp_nulibtable_logitemp_max_;
	nulibtable_logIeta_min			= nulibtable_mp_nulibtable_logieta_min_;
	nulibtable_logIeta_max			= nulibtable_mp_nulibtable_logieta_max_;

#elif defined __GNUC__
	nulibtable_number_species      = __nulibtable_MOD_nulibtable_number_species;
	nulibtable_number_easvariables = __nulibtable_MOD_nulibtable_number_easvariables;
	nulibtable_number_groups       = __nulibtable_MOD_nulibtable_number_groups;
	nulibtable_nrho                = __nulibtable_MOD_nulibtable_nrho;
	nulibtable_ntemp               = __nulibtable_MOD_nulibtable_ntemp;
	nulibtable_nye                 = __nulibtable_MOD_nulibtable_nye;
	nulibtable_nItemp			   = __nulibtable_MOD_nulibtable_nitemp;
	nulibtable_nIeta			   = __nulibtable_MOD_nulibtable_nieta;
	nulibtable_energies            = __nulibtable_MOD_nulibtable_energies;
	nulibtable_ewidths             = __nulibtable_MOD_nulibtable_ewidths;
	nulibtable_ebottom             = __nulibtable_MOD_nulibtable_ebottom;
	nulibtable_etop                = __nulibtable_MOD_nulibtable_etop;
	nulibtable_logrho              = __nulibtable_MOD_nulibtable_logrho;
	nulibtable_logtemp             = __nulibtable_MOD_nulibtable_logtemp;
	nulibtable_ye                  = __nulibtable_MOD_nulibtable_ye;
	nulibtable_logItemp            = __nulibtable_MOD_nulibtable_logitemp;
	nulibtable_logIeta             = __nulibtable_MOD_nulibtable_logieta;
	nulibtable_logtemp_min         = __nulibtable_MOD_nulibtable_logtemp_min;
	nulibtable_logtemp_max         = __nulibtable_MOD_nulibtable_logtemp_max;
	nulibtable_logrho_min          = __nulibtable_MOD_nulibtable_logrho_min;
	nulibtable_logrho_max          = __nulibtable_MOD_nulibtable_logrho_max;
	nulibtable_ye_min              = __nulibtable_MOD_nulibtable_ye_min;
	nulibtable_ye_max              = __nulibtable_MOD_nulibtable_ye_max;
	nulibtable_logItemp_min		   = __nulibtable_MOD_nulibtable_logitemp_min;
	nulibtable_logItemp_max		   = __nulibtable_MOD_nulibtable_logitemp_max;
	nulibtable_logIeta_min		   = __nulibtable_MOD_nulibtable_logieta_min;
	nulibtable_logIeta_max		   = __nulibtable_MOD_nulibtable_logieta_max;


#endif
}
void nulib_set_globals(){
#ifdef __INTEL_COMPILER
	nulib_total_eos_variables       = nulib_mp_total_eos_variables_;
#elif defined __GNUC__
	nulib_total_eos_variables      = __nulib_MOD_total_eos_variables;
#endif
}

/**********************/
/* nulib_get_nspecies */
/**********************/
int nulib_get_nspecies(){
	assert(nulibtable_number_species > 0);
	return nulibtable_number_species;
}

/**************/
/* nulib_init */
/**************/
void nulib_init(string filename){
	int include_Ielectron = 0;//false;
	int include_epannihil_kernels = 0;//false;
	nulibtable_reader_((char*)filename.c_str(), &include_Ielectron, &include_epannihil_kernels, filename.length());
	nulibtable_set_globals();
}

/**********************/
/* get_nut_eas_arrays */
/**********************/
// leave serial - called per grid zone
void nulib_get_eas_arrays(
		double rho,                     // g/cm^3
		double temp,                    // K
		double ye, int nulibID,
		cdf_array& nut_emiss,         // erg/cm^3/s/ster
		vector<double>& nut_absopac,    // cm^-1
		vector<double>& nut_scatopac){  // cm^-1

	assert(rho >= 0);
	assert(temp >= 0);
	assert(ye >= 0);
	assert(ye <= 1.0);
	assert(nulibID >= 0);

	int nvars    = nulibtable_number_easvariables;
	int ngroups  = nulibtable_number_groups;
	assert(nvars > 0);
	assert(ngroups > 0);

	// apparently it's valid to declare array sizes at runtime like this... if it breaks, use malloc
	double eas_energy[nvars][ngroups]; //[0][*] = energy-bin-integrated emissivity (erg/cm^3/s/ster).
	//[1][*] = absorption opacity (1/cm)
	//[2][*] = scattering opacity (1/cm)

	//check sizes match and we are within the table boundaries
	if(nvars != 3){
		cout << "ERROR: nulibtable has different number of eas variables than 3" << endl;
		exit(EXIT_FAILURE);
	}
	if(ngroups != (int)nut_absopac.size()){
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
		for(int j=0; j<ngroups; j++){
			nut_emiss.set_value(j, eas_energy[0][j]);
			nut_absopac [j] =      eas_energy[1][j];
			nut_scatopac[j] =      eas_energy[2][j];
		}
	}
	nut_emiss.N = NaN;
}


/***********************/
/* nulib_get_pure_emis */
/***********************/
void nulib_get_pure_emis(
		double rho,                     // g/cm^3
		double temp,                    // K
		double ye, int nulibID,
		vector<double>& nut_emiss){   // erg/cm^3/s/ster/Hz

	assert(rho >= 0);
	assert(temp >= 0);
	assert(ye >= 0);
	assert(ye <= 1);
	assert(nulibID >= 0);

	vector<double> widths;
	widths.assign(nulibtable_ewidths, nulibtable_ewidths + nulibtable_number_groups);

	int nvars    = nulibtable_number_easvariables;
	int ngroups  = nulibtable_number_groups;
	assert(nvars > 0);
	assert(ngroups > 0);

	// apparently it's valid to declare array sizes at runtime like this... if it breaks, use malloc
	double eas_energy[nvars][ngroups]; //[0][*] = energy-bin-integrated emissivity (erg/cm^3/s/ster).
	//[1][*] = absorption opacity (1/cm)
	//[2][*] = scattering opacity (1/cm)

	//check sizes match and we are within the table boundaries
	if(nvars != 3){
		cout << "ERROR: nulibtable has different number of eas variables than 3" << endl;
		exit(EXIT_FAILURE);
	}
	if(ngroups != (int)nut_emiss.size()){
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
	nu_grid.x.assign(nulibtable_etop,
			nulibtable_etop + nulibtable_number_groups);

	// convert from MeV to frequency using the Planck constant
	for(unsigned i=0; i<nu_grid.size(); i++) nu_grid.x[i] /= pc::h_MeV;
	nu_grid.min = nulibtable_ebottom[0];
}


/*****************************/
/* get the rho, T, Ye arrays */
/*****************************/
void nulib_get_rho_array(vector<double>& array){ // g/cm^3
	array.assign(nulibtable_logrho,
			nulibtable_logrho  + nulibtable_nrho);
	for(unsigned i=0; i<array.size(); i++) array[i] = pow(10.0, array[i]);
}
void nulib_get_T_array(vector<double>& array){ // K
	array.assign(nulibtable_logtemp, nulibtable_logtemp + nulibtable_ntemp);
	for(unsigned i=0; i<array.size(); i++) array[i] = pow(10.0, array[i]) / pc::k_MeV; // convert from MeV to K
}
void nulib_get_Ye_array(vector<double>& array){
	array.assign(nulibtable_ye, nulibtable_ye + nulibtable_nye);
}
void nulib_get_IT_array(vector<double>& array){ //K
	array.assign(nulibtable_logItemp, nulibtable_logItemp + nulibtable_nItemp);
	for(unsigned i=0; i<array.size(); i++) array[i] = pow(10.0,array[i]) / pc::k_MeV;
}
void nulib_get_Ieta_array(vector<double>& array){
	array.assign(nulibtable_logIeta, nulibtable_logIeta + nulibtable_nIeta);
	for(unsigned i=0; i<array.size(); i++) array[i] = pow(10.0,array[i]);
}

/*************************/
/* nulib_get_{*min,*max} */
/*************************/
double nulib_get_Tmin()   {return pow(10,min(nulibtable_logtemp_min,nulibtable_logItemp_min)) / pc::k_MeV;} // K
double nulib_get_Tmax()   {return pow(10,max(nulibtable_logtemp_max,nulibtable_logItemp_max)) / pc::k_MeV;} // K
double nulib_get_rhomin() {return pow(10,nulibtable_logrho_min);} // g/ccm
double nulib_get_rhomax() {return pow(10,nulibtable_logrho_max);} // g/ccm
double nulib_get_Yemin()  {return nulibtable_ye_min;}
double nulib_get_Yemax()  {return nulibtable_ye_max;}
bool nulib_in_range(double rho /* g/ccm */, double temp /* K */, double ye){
	return temp>=nulib_get_Tmin()   and temp<=nulib_get_Tmax() and
		   rho >=nulib_get_rhomin() and rho <=nulib_get_rhomax() and
		   ye  >=nulib_get_Yemin()  and ye  <=nulib_get_Yemax();
}

/*************/
/* eos calls */
/*************/
void nulib_eos_read_table(char* eos_filename){
	read_eos_table_(eos_filename);
	nulib_set_globals();
}

double nulib_eos_munue(const double rho /* g/ccm */, const double temp /* K */, const double ye){ // erg
	assert(rho>=0);
	assert(temp>=0);
	assert(ye>=0 && ye<=1);
	if(!nulib_in_range(rho,temp,ye)) return 0.0;
	else{
		double eos_variables[nulib_total_eos_variables];
		for(int i=0; i<nulib_total_eos_variables; i++) eos_variables[i] = 0;
		eos_variables[0] = rho;
		eos_variables[1] = temp * pc::k_MeV;
		eos_variables[2] = ye;

		set_eos_variables_(eos_variables);
		double mue = eos_variables[10];
		double muhat = eos_variables[13];

		return (mue-muhat)*pc::MeV_to_ergs;
	}
}

double nulib_eos_mue(const double rho /* g/ccm */, const double temp /* K */, const double ye){ // erg
	assert(rho>=0);
	assert(temp>=0);
	assert(ye>=0 && ye<=1);
	if(!nulib_in_range(rho,temp,ye)) return 0.0;
	else{
		double eos_variables[nulib_total_eos_variables];
		for(int i=0; i<nulib_total_eos_variables; i++) eos_variables[i] = 0;
		eos_variables[0] = rho;
		eos_variables[1] = temp * pc::k_MeV;
		eos_variables[2] = ye;

		set_eos_variables_(eos_variables);
		double mue = eos_variables[10];
		return mue*pc::MeV_to_ergs;
	}
}


/****************************************/
/* Pair production/annihilation kernels */
/****************************************/
void nulib_get_epannihil_kernels(const double rho,                    /* g/ccm */
								 const double temp,                   /* K */
								 const double ye, const int nulibID,
								 vector< vector< vector<double> > >& phi_production,	 /* ccm/s */ // phi[legendre_index][this_group][anti-group]
								 vector< vector< vector<double> > >& phi_annihilation){  /* ccm/s */ // phi[legendre_index][this_group][anti-group]

	assert(temp <= pow(10.0,nulibtable_logItemp_max));
	assert(temp >= pow(10.0,nulibtable_logItemp_min));
	assert(nulibID >= 0);

	// fetch the relevant table from nulib. NuLib only accepts doubles.
	double temp_MeV = temp * pc::k_MeV; // MeV
	int lns = nulibID+1;                // fortran array indices start with 1

	// get the chemical potential
	double munue = nulib_eos_munue(rho,temp,ye); // MeV
	double eta = munue/temp_MeV;
	assert(eta >= pow(10.0,nulibtable_logIeta_min));
	assert(eta <= pow(10.0,nulibtable_logIeta_max));

	// apparently it's valid to declare array sizes at runtime like this... if it breaks, use malloc
	int ngroups = nulibtable_number_groups;
	double eas[4][ngroups][ngroups]; //[a][j][i] = kernels for nulib_ID's group i and the anti-group j (ccm/s).
									 // for a, 0->phi_0^p  1->phi_0^a  2->phi_1^p  3->phi_1^a

	// read the kernels from nulib
	int n_legendre_coefficients = 2;
	int n_phis = 2*n_legendre_coefficients;
	nulibtable_epannihil_single_species_range_energy_(&temp_MeV, &eta, &lns, (double*)eas, &ngroups, &ngroups, &n_phis);

	// assign values to phi_production and phi_annihilation
	phi_production.resize  (n_legendre_coefficients);
	phi_annihilation.resize(n_legendre_coefficients);
	for(int l=0; l<n_legendre_coefficients; l++){
		phi_production  [l].resize(ngroups);
		phi_annihilation[l].resize(ngroups);
		for(int j=0; j<ngroups; j++){
			phi_production  [l][j].resize(ngroups);
			phi_annihilation[l][j].resize(ngroups);
			for(int i=0; i<ngroups; i++){
				phi_production  [l][j][i] = eas[2*l  ][j][i];
				phi_annihilation[l][j][i] = eas[2*l+1][j][i];
			}
		}
	}
}
