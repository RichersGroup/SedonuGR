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

#include <mpi.h>
#include <cmath>
#include <string>
#include <cstdlib>
#include "nulib_interface.h"
#include "global_options.h"
#include "H5Cpp.h"
#include "Axis.h"

using namespace std;
namespace pc = physical_constants;


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
		double* eta,   // mu/kT
		int* lns,      // species number
		double* phi,   // phi[legendre-p/a index][this_group][anti-group]
		int* ngroups1,
		int* ngroups2,
		int* n_phis);
void nulibtable_inelastic_single_species_range_energy_(
		double* temp,  // MeV
		double* eta,   // mu/kT
		int* lns,      // species number
		double* phi,   // phi[legendre index][out group][in group]
		int* ngroups1, // ng in
		int* ngroups2, // ng out (should be same as eas_n1)
		int* n_phis);   // number of legendre terms (=2)

void nulibtable_reader_(char*,int*,int*,int*,int);

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
int     read_Ielectron;
int     read_epannihil;
int     read_delta;
int     output_scattering_kernels;

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
extern double* nulibtable_mp_nulibtable_logitemp_;
extern double* nulibtable_mp_nulibtable_logieta_;
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
#elif defined _CRAYC
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

#elif defined _CRAYC
	// do nothing
#elif defined __PGI
	// do nothing
#endif
}
void nulib_set_globals(){
#ifdef __INTEL_COMPILER
	nulib_total_eos_variables       = nulib_mp_total_eos_variables_;
#elif defined __GNUC__
	nulib_total_eos_variables      = __nulib_MOD_total_eos_variables;
#elif defined _CRAYC
	// do nothing
#elif defined __PGI
	// do nothing
#endif
}

/**********************/
/* nulib_get_nspecies */
/**********************/
int nulib_get_nspecies(){
	PRINT_ASSERT(nulibtable_number_species,>,0);
	return nulibtable_number_species;
}

/**************/
/* nulib_init */
/**************/
void nulib_init(string filename, int use_scattering_kernels, int use_annihil_kernels){
	read_Ielectron = 0;
	read_epannihil = 0;
	read_delta = 0;
	output_scattering_kernels = use_scattering_kernels;
    if(use_scattering_kernels == 0) {
    	// you should be using transport opacities with isotropic scattering
    	PRINT_ASSERT(hdf5_dataset_exists(filename.c_str(),"/scattering_delta"),==,false);
    	// if inelastic kernels exist in the NuLib file, you should be using them.
    	PRINT_ASSERT(hdf5_dataset_exists(filename.c_str(),"/inelastic_phi0"),==,false);
    }
    else{
    	if(hdf5_dataset_exists(filename.c_str(),"/scattering_delta")) read_delta = 1;
    	if(hdf5_dataset_exists(filename.c_str(),"/inelastic_phi0"))   read_Ielectron = 1;
    }

    if(use_annihil_kernels==1){
    	assert(hdf5_dataset_exists(filename.c_str(),"/epannihil_phi0"));
    	read_epannihil = 1;
    }

	nulibtable_reader_((char*)filename.c_str(), &read_Ielectron, &read_epannihil, &read_delta, filename.length());
	nulibtable_set_globals();

	// output some facts about the table
	int my_rank=-1;
	MPI_Comm_rank( MPI_COMM_WORLD, &my_rank );
	if(my_rank==0){
		cout << "#   rho range: {" << nulib_get_rhomin() << "," << nulib_get_rhomax() << "} g/ccm" << endl;
		cout << "#   T   range: {" << nulib_get_Tmin()*pc::k_MeV << "," << nulib_get_Tmax()*pc::k_MeV << "} MeV" << endl;
		cout << "#   Ye  range: {" << nulib_get_Yemin() << "," << nulib_get_Yemax() << "}" << endl;
		cout << "#   E   range: {" << nulibtable_ebottom[0] << "," << nulibtable_etop[nulibtable_number_groups-1] << "} MeV" << endl;
		cout << "#   n_rho   = " << nulibtable_nrho << endl;
		cout << "#   n_T     = " << nulibtable_ntemp << endl;
		cout << "#   n_Ye    = " << nulibtable_nye << endl;
		cout << "#   n_E     = " << nulibtable_number_groups << endl;
		cout << "#   n_Ieta  = " << nulibtable_nIeta << endl;
		cout << "#   n_Itemp = " << nulibtable_nItemp << endl;
	}
}


/********************************/
/* Inelastic scattering kernels */
/********************************/
void nulib_get_iscatter_kernels(
		const double rho, // g/ccm
		const double temp, // K
		const double ye,
		const int nulibID,
		vector<double>& nut_scatopac, // 1/cm   opac[group in] Input AND output
		vector< vector<double> >& phi0,         // 2pi h^-3 c^-4 phi0 delta(E^3/3)/deltaE [group in][group out] units 1/cm Output.
		vector< vector<double> >& scattering_delta){  // 3.*phi1/phi0   [group_in][group_out] Input AND output.

	// fetch the relevant table from nulib. NuLib only accepts doubles.
	double temp_MeV = temp * pc::k_MeV; // MeV
	int lns = nulibID+1;                // fortran array indices start with 1

	//PRINT_ASSERT(nulibtable_nIeta,>,0);
	//PRINT_ASSERT(nulibtable_nItemp,>,0);
	PRINT_ASSERT(nulibID,>=,0);

	// get the chemical potential
	double mue = nulib_eos_mue(rho,temp,ye)*pc::ergs_to_MeV; // MeV
	double eta = mue/temp_MeV;
	eta = max(eta,pow(10.0,nulibtable_logIeta_min));
	eta = min(eta,pow(10.0,nulibtable_logIeta_max));
	PRINT_ASSERT(eta,<=,pow(10.0,nulibtable_logIeta_max));

	int n_legendre_coefficients = 2;
	int ngroups = nulibtable_number_groups;
	double phi[n_legendre_coefficients][ngroups][ngroups]; //[a][j][i] = legendre index a, out group i, and in group j (ccm/s)

	// read the kernels from nulib
	if(read_Ielectron>0){
		PRINT_ASSERT(temp_MeV,<=,pow(10.0,nulibtable_logItemp_max));
		PRINT_ASSERT(temp_MeV,>=,pow(10.0,nulibtable_logItemp_min));
		nulibtable_inelastic_single_species_range_energy_(&temp_MeV, &eta, &lns, (double*)phi,
				&ngroups, &ngroups, &n_legendre_coefficients);
	}

	// set the arrays.
	double constants = pow(pc::h,3) * pow(pc::c,4) / (4.*pc::pi);
	for(int igin=0; igin<nulibtable_number_groups; igin++){
		double elastic_opac = nut_scatopac[igin];
		nut_scatopac[igin] = 0;
		for(int igout=0; igout<nulibtable_number_groups; igout++){
			double E1 = nulibtable_ebottom[igout] * pc::MeV_to_ergs;
			double E2 = nulibtable_etop[   igout] * pc::MeV_to_ergs;
			double dE3 = E2*E2*E2 - E1*E1*E1;
			double coeff = (dE3/3.0) / constants;

			double inelastic_phi0=0, inelastic_phi1=0;
			if(read_Ielectron){
				inelastic_phi0 = phi[0][igout][igin] * coeff;
				inelastic_phi1 = phi[1][igout][igin] * coeff;
				PRINT_ASSERT(inelastic_phi1,<=,inelastic_phi0);
			}
			// add elastic scatter opacity to the kernel
			if(igin == igout){
				inelastic_phi0 += elastic_opac;
				if(nulibtable_number_easvariables==4)
					inelastic_phi1 += elastic_opac * scattering_delta[igin][igin] / 3.0; // scattering_delta set in nulib_get_eas_arrays
			}
			phi0[igin][igout] = inelastic_phi0;
			nut_scatopac[igin] += inelastic_phi0;
			scattering_delta[igin][igout] = (inelastic_phi0==0 ? 0 : 3. * inelastic_phi1 / inelastic_phi0);
			PRINT_ASSERT(abs(scattering_delta[igin][igout]),<=,3.0+TINY);
			scattering_delta[igin][igout] = min(3., max(-3., scattering_delta[igin][igout]));
		}
	}
}

/************************/
/* Annihilation kernels */
/************************/
void nulib_get_epannihil_kernels(
		const double rho, // g/ccm
		const double temp, // K
		const double ye,
		const int nulibID,
		vector< vector< vector<double> > >& phi){ // 2pi h^-3 c^-4 phi0 delta(E^3/3)/deltaE [group in][group out] units 1/cm/erg

	// set up phi vector
	phi.resize(2);
	for(size_t i=0; i<2; i++)
		phi[i] = vector< vector<double> >(nulibtable_number_groups, vector<double>(nulibtable_number_groups,0));

	// fetch the relevant table from nulib. NuLib only accepts doubles.
	double temp_MeV = temp * pc::k_MeV; // MeV
	int lns = nulibID+1;                // fortran array indices start with 1

	//PRINT_ASSERT(nulibtable_nIeta,>,0);
	//PRINT_ASSERT(nulibtable_nItemp,>,0);
	PRINT_ASSERT(nulibID,>=,0);

	// get the chemical potential
	double mue = nulib_eos_mue(rho,temp,ye) * pc::ergs_to_MeV; // MeV
	double eta = mue/temp_MeV;
	eta = max(eta,pow(10.0,nulibtable_logIeta_min));
	eta = min(eta,pow(10.0,nulibtable_logIeta_max));
	PRINT_ASSERT(eta,<=,pow(10.0,nulibtable_logIeta_max));

	int n_legendre_coefficients = 2*2;
	int ngroups = nulibtable_number_groups;
	double phi_tmp[n_legendre_coefficients][ngroups][ngroups]; //[a][j][i] = legendre index a, out group i, and in group j (ccm/s)

	// read the kernels from nulib
	if(read_epannihil>0){
		PRINT_ASSERT(temp_MeV,<=,pow(10.0,nulibtable_logItemp_max));
		PRINT_ASSERT(temp_MeV,>=,pow(10.0,nulibtable_logItemp_min));
		nulibtable_epannihil_single_species_range_energy_(&temp_MeV, &eta, &lns, (double*)phi_tmp,
				&ngroups, &ngroups, &n_legendre_coefficients);
	}

	// set the arrays.
	for(int igin=0; igin<nulibtable_number_groups; igin++){
		for(int igout=0; igout<nulibtable_number_groups; igout++){
			phi[0][igin][igout] = phi_tmp[1][igout][igin];
			phi[1][igin][igout] = phi_tmp[3][igout][igin];
			PRINT_ASSERT(phi[0][igin][igout],>=,0);
			PRINT_ASSERT(abs(phi[1][igin][igout]),<=,phi[0][igin][igout]); // mathematically impossible to have phi1/phi0>9
		}
	}
}

/**********************/
/* get_nut_eas_arrays */
/**********************/
// leave serial - called per grid zone
void nulib_get_eas_arrays(
		double rho,                     // g/cm^3
		double temp,                    // K
		double ye, int nulibID,
		vector<double>& nut_absopac,    // cm^-1
		vector<double>& nut_scatopac,   // cm^-1
		vector< vector<double> >& scattering_phi0, // 2pi h^-3 c^-4 phi0 delta(E^3/3)/deltaE [group in][group out] units 1/cm Output.
		vector< vector<double> >& scattering_delta){

	PRINT_ASSERT(rho,>=,0);
	PRINT_ASSERT(temp,>=,0);
	PRINT_ASSERT(ye,>=,0);
	PRINT_ASSERT(ye,<=,1.0);
	PRINT_ASSERT(nulibID,>=,0);

	int nvars    = nulibtable_number_easvariables;
	int ngroups  = nulibtable_number_groups;
	PRINT_ASSERT(nvars,>,0);
	PRINT_ASSERT(ngroups,>,0);

	// apparently it's valid to declare array sizes at runtime like this... if it breaks, use malloc
	double eas_energy[nvars][ngroups]; //[0][*] = energy-bin-integrated emissivity (erg/cm^3/s/ster).
	//[1][*] = absorption opacity (1/cm)
	//[2][*] = scattering opacity (1/cm)

	//check sizes match and we are within the table boundaries
	if(read_delta) PRINT_ASSERT(nvars,==,4);
	else           PRINT_ASSERT(nvars,==,3);
	PRINT_ASSERT(ngroups,==,(int)nut_absopac.size());

	// fetch the relevant table from nulib. NuLib only accepts doubles.
	double temp_MeV = temp * pc::k_MeV; // MeV
	int lns = nulibID+1;                // fortran array indices start with 1

	// keep ye in bounds
	ye = min(nulibtable_ye_max,ye);
	ye = max(nulibtable_ye_min,ye);

	// If the density or temperature are too low, just set everything to zero
	if(log10(rho) < nulibtable_logrho_min || log10(temp_MeV) < nulibtable_logtemp_min)
		for(int j=0; j<ngroups; j++){
			nut_absopac [j] = 0;
			nut_scatopac[j] = 0;
		}

	// Otherwise, fill with the appropriate values
	else{
		// nulib's emis is bin integrated. Summing these values would correspond to our intended nut_emiss value at the bin top.
		// must rebin to get the integrated value to be at the same location as the opacities. (CDF value corresponds to emission rate at or below that energy)
		nulibtable_single_species_range_energy_(&rho, &temp_MeV, &ye, &lns,
				(double*)eas_energy, &ngroups, &nvars);
		for(int j=0; j<ngroups; j++){
			nut_absopac [j] = eas_energy[1][j];
			nut_scatopac[j] = eas_energy[2][j];
			if(nulibtable_number_easvariables==4)
				scattering_delta[j][j] = eas_energy[3][j];
 		}

		// set inelastic kernels if they exist in the table
		if(output_scattering_kernels!=0)
			nulib_get_iscatter_kernels(rho,temp,ye,nulibID,nut_scatopac,scattering_phi0,scattering_delta);
	}
}


/*********************/
/* nulib_get_nu_grid */
/*********************/
// Fill in the locate array with values of the array stored in the fortran module
void nulib_get_nu_grid(Axis& nu_grid){ // Hz
	vector<double> tops(nulibtable_number_groups), mid(nulibtable_number_groups);

	// assign values from the NuLib module to nu_grid
	tops.assign(nulibtable_etop,
			nulibtable_etop + nulibtable_number_groups);
	mid.assign(nulibtable_energies,
			nulibtable_energies + nulibtable_number_groups);

	// convert from MeV to frequency using the Planck constant
	for(size_t i=0; i<tops.size(); i++){
		tops[i] /= pc::h_MeV;
		mid[i] /= pc::h_MeV;
	}
	double minval = nulibtable_ebottom[0] / pc::h_MeV;

	nu_grid = Axis(minval, tops, mid);
}


/*****************************/
/* get the rho, T, Ye arrays */
/*****************************/
void nulib_get_rho_array(vector<double>& array){ // g/cm^3
	array.assign(nulibtable_logrho,
			nulibtable_logrho  + nulibtable_nrho);
	for(size_t i=0; i<array.size(); i++) array[i] = pow(10.0, array[i]);
}
void nulib_get_T_array(vector<double>& array){ // K
	array.assign(nulibtable_logtemp, nulibtable_logtemp + nulibtable_ntemp);
	for(size_t i=0; i<array.size(); i++) array[i] = pow(10.0, array[i]) / pc::k_MeV; // convert from MeV to K
}
void nulib_get_Ye_array(vector<double>& array){
	array.assign(nulibtable_ye, nulibtable_ye + nulibtable_nye);
}
void nulib_get_IT_array(vector<double>& array){ //K
	array.assign(nulibtable_logItemp, nulibtable_logItemp + nulibtable_nItemp);
	for(size_t i=0; i<array.size(); i++) array[i] = pow(10.0,array[i]) / pc::k_MeV;
}
void nulib_get_Ieta_array(vector<double>& array){
	array.assign(nulibtable_logIeta, nulibtable_logIeta + nulibtable_nIeta);
	for(size_t i=0; i<array.size(); i++) array[i] = pow(10.0,array[i]);
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
	PRINT_ASSERT(rho,>=,0);
	PRINT_ASSERT(temp,>=,0);
	PRINT_ASSERT(ye,>=,0);
	PRINT_ASSERT(ye,<=,1);
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
	PRINT_ASSERT(rho,>=,0);
	PRINT_ASSERT(temp,>=,0);
	PRINT_ASSERT(ye,>=,0);
	PRINT_ASSERT(ye,<=,1);
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


