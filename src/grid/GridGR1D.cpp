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
#include <fstream>
#include "Transport.h"
#include "GridGR1D.h"
#include "global_options.h"

using namespace std;
namespace pc = physical_constants;

GridGR1D::GridGR1D(){
	PRINT_ASSERT(NDIMS,==,1);
	grid_type = "GridGR1D";
	ghosts1=-1;
}

//------------------------------------------------------------
// initialize the zone geometry from model file
//------------------------------------------------------------
void GridGR1D::read_model_file(Lua* lua)
{
	// DO NOTHING - zones will be set by a different function during each iteration.
}

void GridGR1D::init(Lua* lua){
	// DO NOTHING
}

//------------------------------------------------------------
// Reflect off symmetry axis
//------------------------------------------------------------
void GridGR1D::symmetry_boundaries(EinsteinHelper *eh) const{
	// NONE - just flow out of outer boundary
}

void GridGR1D::set_fluid(const double* rho_in, const double* T_in, const double* Ye_in, const double* vr_in){
	for(int z_ind=0; z_ind<rho.size(); z_ind++)
	{
		rho[z_ind]  = rho_in[z_ind+ghosts1];
		T[z_ind]    =   T_in[z_ind+ghosts1];
		Ye[z_ind]   =  Ye_in[z_ind+ghosts1];
		vr[z_ind]   =  vr_in[z_ind+ghosts1];

		H_vis[z_ind] = 0;
		PRINT_ASSERT(xAxes[0].top[z_ind],>,(z_ind==0 ? xAxes[0].min : xAxes[0].top[z_ind-1]));
		PRINT_ASSERT(rho[z_ind],>=,0);
		PRINT_ASSERT(T[z_ind],>=,0);
		PRINT_ASSERT(Ye[z_ind],>=,0);
		PRINT_ASSERT(Ye[z_ind],<=,1.0);
	}
}

void GridGR1D::initialize_grid(const double* rarray, const int n_zones, const int nghost){
	ghosts1 = nghost;
	vector<double> rtop(n_zones), rmid(n_zones);
	double rmin=0;
	for(int z_ind=0; z_ind<n_zones; z_ind++){
		rtop[z_ind] = rarray[z_ind + ghosts1+1];
		double last = z_ind==0 ? rmin : rtop[z_ind-1];
		PRINT_ASSERT(rtop[z_ind],>,last);
		rmid[z_ind] = 0.5 * (rtop[z_ind] + last);
	}
	xAxes[0] = Axis(rmin, rtop, rmid);

	rho.set_axes(xAxes);
	T.set_axes(xAxes);
	Ye.set_axes(xAxes);
	H_vis.set_axes(xAxes);

	cout << "#   Sedonu grid has "<< rho.size() << " zones." << endl;
	cout << "#   Sedonu outer boundary is at "<< xAxes[0].top[n_zones-1] << " cm" << endl;
}
