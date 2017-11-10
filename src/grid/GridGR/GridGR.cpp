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
	lh->particle_copy(lab).print();

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

	// make vector null again
	normalize_null(knew,4,xnew);

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
double GridGR::dot3(const double a[3], const double b[3], const int size, const double xup[4]) const{
	PRINT_ASSERT(size,==,3);

	double product = 0;
	for(int mu=0; mu<3; mu++) for(int nu=0; nu<3; nu++)
		product += a[mu] * b[nu] * g_down(xup,mu,nu);
	return product;
}
void GridGR::normalize(double a[4], const int size, const double xup[4]) const{
	PRINT_ASSERT(size,==,4);

	double inv_norm = 0;
	for(int mu=0; mu<4; mu++) for(int nu=0; nu<4; nu++)
		inv_norm += a[mu] * a[nu] * g_down(xup,mu,nu);
	inv_norm = 1.0 / sqrt(inv_norm);

	for(int mu=0; mu<4; mu++) a[mu] *= inv_norm;
}
void GridGR::normalize_null(double a[4], const int size, const double xup[4]) const{
	PRINT_ASSERT(size,==,4);

	double A = g_down(xup,3,3);

	double B = 0;
	for(int i=0; i<3; i++) B += a[i] * g_down(xup,i,3);
	B *= 2.0;

	double C = 0;
	for(int i=0; i<3; i++) for(int j=0; j<3; j++) C += a[i] * a[j] * g_down(xup,i,j);

	a[3] = (-B - sqrt(abs( B*B - 4.0*A*C )) ) / (2.0*A);
}

void GridGR::orthogonalize(double v[4], const double e[4], const double xup[4], const int size) const{
	PRINT_ASSERT(abs(dot(e,e,size,xup)-1.0),<,TINY); // assume the basis vector is normalized
	PRINT_ASSERT(size,==,4);

	double projection = dot(v,e,4,xup);
	for(int mu=0; mu<4; mu++) v[mu] -= projection * v[mu];
}
void GridGR::tetrad_to_coord(const double xup[4], const double u[4], double kup_tetrad[4], const int size) const{
	PRINT_ASSERT(size,==,4);
	double e[4][4];

	// normalize four-velocity to get timelike vector
	for(int mu=0; mu<4; mu++) e[3][mu] = u[mu];
	normalize(e[3], 4, xup);

	// use x0 as a trial vector
	e[0][0] = 1.0;
	e[0][1] = 0;
	e[0][2] = 0;
	e[0][3] = 0;
	orthogonalize(e[0],e[3],xup,4);
	normalize(e[0],4,xup);

	// use x1 as a trial vector
	e[1][0] = 0;
	e[1][1] = 1.0;
	e[1][2] = 0;
	e[1][3] = 0;
	orthogonalize(e[1],e[3],xup,4);
	orthogonalize(e[1],e[0],xup,4);
	normalize(e[1],4,xup);

	// use x2 as a trial vector
	e[2][0] = 0;
	e[2][1] = 0;
	e[2][2] = 1.0;
	e[2][3] = 0;
	orthogonalize(e[2],e[3],xup,4);
	orthogonalize(e[2],e[0],xup,4);
	orthogonalize(e[2],e[1],xup,4);
	normalize(e[2],4,xup);

	// sanity checks
	PRINT_ASSERT(abs(dot(e[0],e[1],4,xup)),<,TINY);
	PRINT_ASSERT(abs(dot(e[0],e[2],4,xup)),<,TINY);
	PRINT_ASSERT(abs(dot(e[0],e[3],4,xup)),<,TINY);
	PRINT_ASSERT(abs(dot(e[1],e[2],4,xup)),<,TINY);
	PRINT_ASSERT(abs(dot(e[1],e[3],4,xup)),<,TINY);
	PRINT_ASSERT(abs(dot(e[2],e[3],4,xup)),<,TINY);

	// transform to coordinate frame
	double kup[4];
	for(int mu=0; mu<4; mu++) kup[mu] = dot(kup_tetrad,e[mu],4,xup);
	for(int mu=0; mu<4; mu++) kup_tetrad[mu] = kup[mu];
}

// isotropic scatter, done in COMOVING frame
void GridGR::isotropic_kup(const double nu, double kup[4], const double xup[4], const int size, ThreadRNG *rangen) const
{
	PRINT_ASSERT(size,==,4);

	// Randomly generate new direction isotropically in comoving frame
	double D[3];
	isotropic_direction(D,3,rangen);

	kup[0] = nu * D[0] * 2.0*pc::pi / pc::c;
	kup[1] = nu * D[1] * 2.0*pc::pi / pc::c;
	kup[2] = nu * D[2] * 2.0*pc::pi / pc::c;
	kup[3] = nu        * 2.0*pc::pi / pc::c;

	// move from tetrad frame to comoving frame
	double v[3];
	interpolate_fluid_velocity(xup,4,v,3);
	double u[4];
	double gamma = LorentzHelper::lorentz_factor(v,3);
	u[0] = gamma * v[0];
	u[1] = gamma * v[1];
	u[2] = gamma * v[2];
	u[3] = gamma * pc::c;
	tetrad_to_coord(xup, u, kup, 4);
}
