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

#ifndef _NULIB_INTERFACE_H
#define _NULIB_INTERFACE_H

#include <vector>
#include "CDFArray.h"
#include "Axis.h"

using namespace std;
//
// returns everything in standard CGS units (i.e. ergs, s, cm, K, Hz)

void nulib_init(string filename, int include_Ielectron, int do_annihil);
void nulib_get_eas_arrays(double rho, double temp, double ye, int nulibID,
		vector<double>& nut_BB, vector<double>& nut_absopac, vector<double>& nut_scatopac,
		vector<CDFArray>& EoutCDF, vector< vector<double> >& phi1_phi0);
void nulib_get_epannihil_kernels(
		const double rho, const double temp, const double ye, const int nulibID,
		vector< vector< vector<double> > >& phi);
void nulib_get_pure_emis(double rho, double temp, double ye, int nulibID, vector<double>& pure_emis);
void nulib_get_nu_grid(Axis& nut_nu_grid);
int nulib_get_nspecies();

double nulib_get_Tmin();
double nulib_get_Tmax();
double nulib_get_Yemin();
double nulib_get_Yemax();
double nulib_get_rhomin();
double nulib_get_rhomax();
bool nulib_in_range(double rho /* g/ccm */, double temp /* K */, double ye);

void   nulib_eos_read_table(char* eos_filename);
double nulib_eos_munue(const double rho, const double temp, const double ye);
double nulib_eos_mue(const double rho, const double temp, const double ye);

void nulib_get_rho_array(vector<double>& array);
void nulib_get_T_array(vector<double>& array);
void nulib_get_Ye_array(vector<double>& array);
#endif
