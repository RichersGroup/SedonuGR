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

#ifndef _NEUTRINOS_H
#define _NEUTRINOS_H

#include "species_general.h"

class neutrinos: public species_general
{

protected:

public:

	neutrinos();
	virtual ~neutrinos() {}

	int num_nut_species;
	int nulibID;
	double cutoff;

	// required functions
	void myInit(Lua* lua);
	void set_eas(int zone_index);
	double blackbody(const double T, const double chempot, const double nu) const;
	static double annihilation_rate(const spectrum_array& nu_dist, const spectrum_array& nbar_dist, const bool electron_type, const int weight);
};

#endif
