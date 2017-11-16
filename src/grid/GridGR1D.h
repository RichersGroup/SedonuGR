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

#ifndef _GRID_GR1D_H
#define _GRID_GR1D_H 1

#include "Grid1DSphere.h"

//*******************************************
// 1-Dimensional Spherical geometry
//*******************************************
class GridGR1D: public Grid1DSphere
{

private:
	// store location of the outer edge of the zone.
	//LocateArray r_out;
	int ghosts1;

public:

	GridGR1D();
	virtual ~GridGR1D() {}


	// required functions
	void read_model_file(Lua* lua);
	void symmetry_boundaries      (LorentzHelper *lh) const;
	void init(Lua* lua);

	// GR1D-specific functions
	void initialize_grid(const double* x, const int length, const int nghost);
	void set_fluid(const double* rho, const double* T, const double* Ye, const double* vr);
};


#endif
