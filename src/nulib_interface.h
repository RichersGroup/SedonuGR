#ifndef _NULIB_INTERFACE_H
#define _NULIB_INTERFACE_H
#include <vector>
#include "cdf_array.h"
#include "locate_array.h"

using namespace std;

void nulib_init(string filename);
void nulib_get_eas_arrays(real rho, real temp, real ye, int nulibID,
			      cdf_array& nut_emiss, vector<real>& nut_absopac, vector<real>& nut_scatopac);
void nulib_get_nu_grid(locate_array& nut_nu_grid);
int nulib_get_nspecies();

double nulib_get_Tmin();
double nulib_get_Tmax();
double nulib_get_Yemin();
double nulib_get_Yemax();

#endif
