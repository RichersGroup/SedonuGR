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

#ifndef _GRID_GR_H
#define _GRID_GR_H 1

#include "Grid.h"

using namespace std;

enum Index{up,down};

//*******************************************
// general GR grid functionality
//*******************************************
class GridGR: public Grid
{

private:

	class Metric{
	private:
		double gup[4][4];
		double gdown[4][4];
	public:
		template<Index updown>
		double dot3(const double a[4], const double b[4], const int size) const;

		template<Index updown>
		double dot4(const double a[4], const double b[4], const int size) const;
	};

	class Christoffel{
	private:
		double gamma[4][4][4];
	public:
		double contract(const int upindex, const double p1up[4], const double p2up[4], const int size);
	};

	vector<Metric> metric;

public:

	virtual ~GridGR() {}

};


#endif
