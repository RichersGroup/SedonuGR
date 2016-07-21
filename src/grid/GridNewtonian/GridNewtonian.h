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

#ifndef _GRID_NEWT_H
#define _GRID_NEWT_H 1

#include "Grid.h"

using namespace std;

//*******************************************
// general Newtonian grid functionality
//*******************************************
class GridNewtonian: public Grid
{

private:

public:

	virtual ~GridNewtonian() {}


	void integrate_geodesic(LorentzHelper *lh) const;
	void random_core_x_D(const double r_core, ThreadRNG *rangen, double x3[3], double D[3], const int size) const;

	// vector functions
	double dot(const double a[], const double b[], const int size, const double xup[]) const;
	void normalize(double a[], const int size, const double xup[]) const;
	void normalize_null(double a[], const int size, const double xup[]) const;
};


#endif
