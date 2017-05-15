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
	grid_type = "GridGR1D";
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
void GridGR1D::symmetry_boundaries(LorentzHelper *lh) const{
	// NONE - just flow out of outer boundary
}

void GridGR1D::set_fluid(const double* rho, const double* T, const double* Ye, const double* vr, const double* vphi, const int n_zones){
	PRINT_ASSERT(n_zones,==,z.size());
	for(int z_ind=0; z_ind<n_zones; z_ind++)
	{
		z[z_ind].rho = rho[z_ind];
		z[z_ind].T = T[z_ind];
		z[z_ind].Ye = Ye[z_ind];
		z[z_ind].u[0] = vr[z_ind];
		z[z_ind].u[1] = vphi[z_ind];
		z[z_ind].u[2] = 0;

		z[z_ind].H_vis = 0;
		PRINT_ASSERT(r_out[z_ind],>,(z_ind==0 ? r_out.min : r_out[z_ind-1]));
		PRINT_ASSERT(z[z_ind].rho,>=,0);
		PRINT_ASSERT(z[z_ind].T,>=,0);
		PRINT_ASSERT(z[z_ind].Ye,>=,0);
		PRINT_ASSERT(z[z_ind].Ye,<=,1.0);
		PRINT_ASSERT(z[z_ind].u[0]*z[z_ind].u[0] + z[z_ind].u[1]*z[z_ind].u[1] + z[z_ind].u[2]*z[z_ind].u[2],<,pc::c*pc::c);
	}
}

void GridGR1D::initialize_grid(const double* rarray, const int n_zones){
	z.resize(n_zones);
	r_out.resize(n_zones);
	r_out.min = 0;
	for(int z_ind=0; z_ind<n_zones; z_ind++){
		r_out[z_ind] = rarray[z_ind];
		if(z_ind==0) PRINT_ASSERT(r_out[z_ind],>,0);
		else PRINT_ASSERT(r_out[z_ind],>,r_out[z_ind-1]);
	}

	cout << "#   Sedonu grid has "<< z.size() << " zones." << endl;
	cout << "#   Sedonu outer boundary is at "<< r_out[n_zones-1] << " km" << endl;
}
