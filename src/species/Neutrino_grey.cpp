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
#include "Neutrino_grey.h"
#include "Transport.h"
#include "Grid.h"
#include "nulib_interface.h"

using namespace std;
namespace pc = physical_constants;


Neutrino_grey::Neutrino_grey(){
	grey_opac = NaN;
	grey_abs_frac = NaN;
}

//----------------------------------------------------------------
// called from species_general::init (neutrino-specific stuff)
//----------------------------------------------------------------
void Neutrino_grey::myInit(Lua* lua)
{
	grey_abs_frac = lua->scalar<double>("grey_abs_frac");
	grey_opac     = lua->scalar<double>("grey_opac");

	// let the base class do the rest
	Neutrino::myInit(lua);
}


//-----------------------------------------------------------------
// set emissivity, abs. opacity, and scat. opacity in zones
//-----------------------------------------------------------------
void Neutrino_grey::set_eas(int zone_index)
{
	Zone* z = &(sim->grid->z[zone_index]);
	double ngroups = (double)emis[zone_index].size();

	PRINT_ASSERT(grey_abs_frac,>=,0);
	PRINT_ASSERT(grey_abs_frac,<=,1.0);
	for(unsigned j=0;j<nu_grid.size();j++)
	{
		double nu  = nu_grid.center(j);        // (Hz)
		double dnu = nu_grid.delta(j);         // (Hz)
		double bb  = blackbody(z->T,0*pc::MeV_to_ergs,nu)*dnu;  // (erg/s/cm^2/ster)

		double a = grey_opac*z->rho*grey_abs_frac;
		double s = grey_opac*z->rho*(1.0-grey_abs_frac);

		emis[zone_index].set_value(j,a*bb); // (erg/s/cm^3/ster)
		abs_opac[zone_index][j] = a;        // (1/cm)
		scat_opac[zone_index][j] = s;       // (1/cm)

		// set the biased emissivity
		biased_emis[zone_index].set_value(j, emis[zone_index].get_value(j)
				* sim->importance(a, s, sim->grid->zone_min_length(zone_index)));
	}

	emis[zone_index].normalize(cutoff/(double)ngroups);
	biased_emis[zone_index].normalize(cutoff/(double)ngroups);
}

//-----------------------------------------------------------------
// get opacity at the frequency
//-----------------------------------------------------------------
void Neutrino_grey::get_opacity(const double com_nu, const int z_ind, double* a, double* s) const
{
	PRINT_ASSERT(z_ind,>=,-1);
	PRINT_ASSERT(com_nu,>,0);

	*a = grey_opac * grey_abs_frac;
	*s = grey_opac * (1.0 - grey_abs_frac);

	PRINT_ASSERT(*a,>=,0);
	PRINT_ASSERT(*s,>=,0);
	PRINT_ASSERT(*a,<,INFINITY);
	PRINT_ASSERT(*s,<,INFINITY);
}
