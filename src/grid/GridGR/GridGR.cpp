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

void GridGR::integrate_geodesic(LorentzHelper *lh){
	double lambda = lh->distance(lab);

	double dk_dlambda[4];
	for(int a=0; a<4; a++){
		dk_dlambda[a] = 0;
		for(int mu=0; mu<4; mu++) for(int nu=0; nu<4; nu++){
			dk_dlambda[a] -= connection_coefficient(lh->p_xup(),a,mu,nu) * lh->p_kup(lab)[mu] * lh->p_kup(lab)[nu];
		}
	}

	// move x forward
	double xnew[4];
	for(int i=0; i<4; i++) xnew[i] = lh->p_xup()[i] + lh->p_kup(lab)[i]*lambda + 0.5*dk_dlambda[i]*lambda*lambda;
	lh->set_p_xup(xnew,4);

	// change k
	//double knew[4];
	//for(int i=0; i<4; i++) knew[i] = lh->p_kup(lab)[i] + dk_dlambda[i]*lambda;
	//lh->set_p_kup<lab>(knew,4);
}
