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

#include "global_options.h"
#include "species_general.h"
#include "transport.h"
#include "locate_array.h"

//-----------------------------------------------------------------
// get opacity at the frequency
//-----------------------------------------------------------------
void species_general::get_opacity(const double com_nu, const int z_ind, double* opac, double* abs_frac) const
{
	PRINT_ASSERT(z_ind,>=,-1);
	PRINT_ASSERT(com_nu,>,0);

	if(z_ind == -1){ // particle is within inner boundary
		*opac = 0;
		*abs_frac = 0;
	}

	else{
		// absorption and scattering opacities
		double a = max(nu_grid.value_at(com_nu, abs_opac[z_ind]),0.0);
		double s = max(nu_grid.value_at(com_nu,scat_opac[z_ind]),0.0);

		// output - net opacity
		*opac = a+s;
		PRINT_ASSERT(*opac,>=,0);

		// output - absorption fraction
		*abs_frac = (a+s>0 ? a/(a+s) : 0);
		PRINT_ASSERT(0,<=,*abs_frac);
		PRINT_ASSERT(*abs_frac,<=,1.0);
	}
}

double species_general::sum_opacity(const int z_ind, const int group) const{
	return abs_opac[z_ind][group] + scat_opac[z_ind][group];
}
