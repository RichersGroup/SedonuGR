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

#ifndef _ZONE_H
#define _ZONE_H

#include <vector>
#include <memory>
#include "global_options.h"

//-------------------------------------------------
// Class to store properties of one zone
//-------------------------------------------------
class Zone
{

public:

	// constructors
	Zone(){
		rho = NaN;
		T = NaN;
		Ye = NaN;
		H_vis = NaN;
		e_abs = NaN;
		nue_abs = NaN;
		anue_abs = NaN;
		e_emit = NaN;
		l_emit = NaN;
	}

	// fluid properties (rho,T are in comoving frame. Ye is invariant.)
	double rho;             // density (g/cm^3)
	double T;               // gas temperature (K)
	double Ye;              // electron fraction
	double H_vis;           // specific heating rate (erg/s/g)

	// radiation quantities (all in comoving frame) (dVdt and lepton number are relativistic invariants)
	double e_abs;                         // radiation energy deposition density rate (ergs/cm^3/s) (comoving frame)
	double nue_abs;                       // electron neutrino number deposition density rate (cm^-3 s^-1) (comoving frame)
	double anue_abs;                      // electron anti-neutrino number deposition density rate (cm^-3 s^-1) (comoving frame)
	double e_emit;                        // radiation energy emission rate (erg/ccm/s) (comoving frame)
	double l_emit;						  // lepton number emission rate (cm^-3 s^-1) (comoving frame)
	std::vector<double> Edens_com;
	std::vector<double> Ndens_com;
};

#endif

// NOTES
// - tried implementing as struct of vector rather than vector of structs (which would eliminate the need for the mpi datatype definitions). Significantly slower, presumably b/c access patterns are not sequential. Also, this probably increased false sharing, making OpenMP less efficient (guess)
