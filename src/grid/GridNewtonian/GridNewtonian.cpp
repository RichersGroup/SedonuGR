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
#include "Lua.h"
#include "nulib_interface.h"
#include "global_options.h"
#include "Transport.h"
#include "H5Cpp.h"
#include "GridNewtonian.h"

using namespace std;
namespace pc = physical_constants;

void GridNewtonian::integrate_geodesic(LorentzHelper *lh) const{
	double xnew[4];
	double Dlab[3];
	lh->p_D(lab,Dlab,3);
	for(int i=0; i<3; i++) xnew[i] = lh->p_xup()[i] + lh->distance(lab) * Dlab[i];
	xnew[3] = lh->p_xup()[3] + lh->distance(lab);
	lh->set_p_xup(xnew,4);
}

// get a random position on the surface of the core
void GridNewtonian::random_core_x_D(const double r_core, ThreadRNG *rangen, double x3[3], double D[3], const int size) const{
	PRINT_ASSERT(size,==,3);

	double phi_core   = 2*pc::pi*rangen->uniform();
	double cosp_core  = cos(phi_core);
	double sinp_core  = sin(phi_core);
	double cost_core  = 1 - 2.0*rangen->uniform();
	double sint_core  = sqrt(1-cost_core*cost_core);

	double a_phot = r_core + r_core*1e-10;
	x3[0] = a_phot*sint_core*cosp_core;
	x3[1] = a_phot*sint_core*sinp_core;
	x3[2] = a_phot*cost_core;

	// pick photon propagation direction wtr to local normal
	double phi_loc = 2*pc::pi*rangen->uniform();
	// choose sqrt(R) to get outward, cos(theta) emission
	double cost_loc  = sqrt(rangen->uniform());
	double sint_loc  = sqrt(1 - cost_loc*cost_loc);
	// local direction vector
	double D_xl = sint_loc*cos(phi_loc);
	double D_yl = sint_loc*sin(phi_loc);
	double D_zl = cost_loc;
	// apply rotation matrix to convert D vector into overall frame
	D[0] = cost_core*cosp_core*D_xl-sinp_core*D_yl+sint_core*cosp_core*D_zl;
	D[1] = cost_core*sinp_core*D_xl+cosp_core*D_yl+sint_core*sinp_core*D_zl;
	D[2] = -sint_core*D_xl+cost_core*D_zl;
	Grid::normalize(D,3);
}
