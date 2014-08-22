#ifndef _NULIB_INTERFACE_H
#define _NULIB_INTERFACE_H
#include <vector>
#include "cdf_array.h"
#include "locate_array.h"
#include "global_options.h"

// returns everything in standard CGS units (i.e. ergs, s, cm, K, Hz)

using namespace std;

void nulib_init(string filename);
void nulib_get_eas_arrays(double rho, double temp, double ye, int nulibID,
		cdf_array& nut_emiss, vector<double>& nut_absopac, vector<double>& nut_scatopac);
void nulib_get_pure_emis(double rho, double temp, double ye, int nulibID, vector<double>& pure_emis);
void nulib_get_nu_grid(locate_array& nut_nu_grid);
int nulib_get_nspecies();

double nulib_get_Tmin();
double nulib_get_Tmax();
double nulib_get_Yemin();
double nulib_get_Yemax();
double nulib_get_rhomin();
double nulib_get_rhomax();

void nulib_get_rho_array(vector<double>& array);
void nulib_get_T_array(vector<double>& array);
void nulib_get_Ye_array(vector<double>& array);
#endif
