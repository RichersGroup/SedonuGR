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
#include "Neutrino_NuLib.h"
#include "Transport.h"
#include "Grid.h"
#include "nulib_interface.h"

using namespace std;
namespace pc = physical_constants;


// constructor
Neutrino_NuLib::Neutrino_NuLib(){
}

//----------------------------------------------------------------
// called from species_general::init (neutrino-specific stuff)
//----------------------------------------------------------------
void Neutrino_NuLib::myInit(Lua* lua)
{
	// set up the frequency table
	nulib_get_nu_grid(nu_grid);

	// set neutrino's min and max values
	T_min  =  nulib_get_Tmin();
	T_max  =  nulib_get_Tmax();
	Ye_min =  nulib_get_Yemin();
	Ye_max =  nulib_get_Yemax();
	rho_min = nulib_get_rhomin();
	rho_max = nulib_get_rhomax();

	// let the base class do the rest
	Neutrino::myInit(lua);
}


//-----------------------------------------------------------------
// set emissivity, abs. opacity, and scat. opacity in zones
//-----------------------------------------------------------------
void Neutrino_NuLib::set_eas(int zone_index)
{
	Zone* z = &(sim->grid->z[zone_index]);
	double ngroups = (double)emis[zone_index].size();

	nulib_get_eas_arrays(z->rho, z->T, z->Ye, ID, cutoff,
			emis[zone_index], abs_opac[zone_index], scat_opac[zone_index]);

	// set the biased emissivity
	for(int g=0; g<(int)nu_grid.size(); g++)
		biased_emis[zone_index].set_value(g, emis[zone_index].get_value(g)
				* sim->importance(abs_opac[zone_index][g], scat_opac[zone_index][g], sim->grid->zone_min_length(zone_index)));

	emis[zone_index].normalize(cutoff/(double)ngroups);
	biased_emis[zone_index].normalize(cutoff/(double)ngroups);
}
