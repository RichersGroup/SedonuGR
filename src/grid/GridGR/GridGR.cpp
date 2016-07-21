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
#include <iostream>
#include <iomanip>
#include <fstream>
#include "GridGR.h"
#include "Lua.h"
#include "nulib_interface.h"
#include "global_options.h"
#include "Transport.h"
#include "H5Cpp.h"

using namespace std;
namespace pc = physical_constants;

void GridGR::integrate_geodesic(LorentzHelper *lh) const{
	double lambda = lh->distance(lab) * pc::c / (2.0*pc::pi * lh->p_nu(lab));

	double dk_dlambda[4];
	for(int a=0; a<4; a++){
		dk_dlambda[a] = 0;
		for(int mu=0; mu<4; mu++) for(int nu=0; nu<4; nu++){
			dk_dlambda[a] -= connection_coefficient(lh->p_xup(),a,mu,nu) * lh->p_kup(lab)[mu] * lh->p_kup(lab)[nu];
		}
		dk_dlambda[a] *= 2.0*pc::pi * lh->p_nu(lab) / pc::c;
	}

	// get new k
	double knew[4];
	for(int i=0; i<4; i++) knew[i] = lh->p_kup(lab)[i] + dk_dlambda[i]*lambda;

	// get new x
	double xnew[4];
	for(int i=0; i<4; i++) xnew[i] = lh->p_xup()[i] + lambda * 0.5 * (lh->p_kup(lab)[i] + knew[i]);

	// actually change the values
	lh->set_p_xup(xnew,4);
	lh->set_p_kup<lab>(knew,4);
}

double GridGR::dot(const double a[4], const double b[4], const int size, const double xup[4]) const{
	PRINT_ASSERT(size,==,4);

	double product = 0;
	for(int mu=0; mu<4; mu++) for(int nu=0; nu<4; nu++)
		product += a[mu] * b[nu] * g_down(xup,mu,nu);
	return product;
}
void GridGR::normalize(double a[], const int size, const double xup[]) const{
	PRINT_ASSERT(size,==,4);

	double inv_norm = 0;
	for(int mu=0; mu<4; mu++) for(int nu=0; nu<4; nu++)
		inv_norm += a[mu] * a[nu] * g_down(xup,mu,nu);
	inv_norm = 1.0 / sqrt(inv_norm);

	for(int mu=0; mu<4; mu++) a[mu] *= inv_norm;
}
void GridGR::normalize_null(double a[], const int size, const double xup[]) const{
	PRINT_ASSERT(size,==,4);

	double A = g_down(xup,3,3);

	double B = 0;
	for(int i=0; i<3; i++) B += a[i] * g_down(xup,i,3);
	B *= 2.0;

	double C = 0;
	for(int i=0; i<3; i++) for(int j=0; j<3; j++) C += a[i] * a[j] * g_down(xup,i,j);

	a[0] = (-B - sqrt(abs( B*B - 4.0*A*C )) ) / (2.0*A);
}
