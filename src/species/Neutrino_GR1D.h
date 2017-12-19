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

#ifndef _NEUTRINOS_GR1D_H
#define _NEUTRINOS_GR1D_H

#include "Neutrino.h"

class Neutrino_GR1D: public Neutrino
{

protected:

	const double nulib_emissivity_gf = 5.59424238e-55/(6.77140812e-6 * 6.77140812e-06 * 6.77140812e-06 * 2.03001708e5);
	const double nulib_opacity_gf = 1.0/6.77140812e-6;
	const double nulib_energy_gf = 1.60217733e-6*5.59424238e-55;
	const double nulib_kernel_gf = 6.77140812e-6*6.77140812e-6*6.77140812e-6/2.03001708e5;

public:

	int ghosts1;
	int n_GR1D_zones;
	int GR1D_emit_outside_shock;
	Neutrino_GR1D();

	void myInit(Lua* lua);
	void set_eas(int zone_index);
	void set_eas_external(const double* easarray, const double GR1D_tau_crit, bool* extract_MC, const double rshock);
	static void set_nu_grid(Lua* lua, LocateArray* nu_grid);
};

#endif
