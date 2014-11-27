#ifndef _NULIB_INTERFACE_H
#define _NULIB_INTERFACE_H
#include "global_options.h"
#include <vector>
#include "cdf_array.h"
#include "locate_array.h"

// returns everything in standard CGS units (i.e. ergs, s, cm, K, Hz)

using namespace std;

void nulib_init(string filename);
void nulib_get_eas_arrays(double rho, double temp, double ye, int nulibID,
		cdf_array& nut_emiss, vector<double>& nut_absopac, vector<double>& nut_scatopac);
void nulib_get_epannihil_kernels(double rho, double temp, double ye, int nulibID,
								 vector< vector< vector<double> > >& phi_production,
								 vector< vector< vector<double> > >& phi_annihilation);
void nulib_get_pure_emis(double rho, double temp, double ye, int nulibID, vector<double>& pure_emis);
void nulib_get_nu_grid(locate_array& nut_nu_grid);
int nulib_get_nspecies();

double nulib_get_Tmin();
double nulib_get_Tmax();
double nulib_get_Yemin();
double nulib_get_Yemax();
double nulib_get_rhomin();
double nulib_get_rhomax();

void   nulib_eos_read_table(char* eos_filename);
double nulib_eos_munue(const double rho, const double temp, const double ye);

void nulib_get_rho_array(vector<double>& array);
void nulib_get_T_array(vector<double>& array);
void nulib_get_Ye_array(vector<double>& array);
#endif
