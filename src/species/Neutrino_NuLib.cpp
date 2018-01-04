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
	// let the base class do the rest
	Neutrino::myInit(lua);
}


//-----------------------------------------------------------------
// set emissivity, abs. opacity, and scat. opacity in zones
//-----------------------------------------------------------------
void Neutrino_NuLib::set_eas(int z_ind)
{
	Zone* z = &(sim->grid->z[z_ind]);
	double ngroups = (double)emis[z_ind].size();

	nulib_get_eas_arrays(sim->grid->rho[z_ind], sim->grid->T[z_ind], sim->grid->Ye[z_ind], ID,
			emis[z_ind], abs_opac[z_ind],
			scat_opac[z_ind], normalized_phi0[z_ind], scattering_delta[z_ind]);

	emis[z_ind].normalize();
}
