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
	unsigned ngroups = sim->grid->nu_grid_axis.size();
	unsigned dir_ind[NDIMS+2];
	sim->grid->rho.indices(z_ind,dir_ind);

	vector<double> tmp_absopac(ngroups), tmp_scatopac(ngroups), tmp_BB(ngroups);
	vector< vector<double> > tmp_delta(ngroups, vector<double>(ngroups));
	vector< vector<double> > tmp_phi0(ngroups, vector<double>(ngroups));
	nulib_get_eas_arrays(sim->grid->rho[z_ind], sim->grid->T[z_ind], sim->grid->Ye[z_ind], ID,
			tmp_BB, tmp_absopac, tmp_scatopac, tmp_phi0, tmp_delta);

	for(unsigned ig=0; ig<ngroups; ig++){
		dir_ind[NDIMS] = ig;
		unsigned global_index = sim->grid->abs_opac[ID].direct_index(dir_ind);
		sim->grid->abs_opac[ID][global_index] = tmp_absopac[ig];
		sim->grid->scat_opac[ID][global_index] = tmp_scatopac[ig];
		sim->grid->BB[ID][global_index] = tmp_BB[ig]  /(pc::h*sim->grid->nu_grid_axis.mid[ig]) * pc::c*pc::c/(4.*pc::pi * sim->grid->nu_grid_axis.delta3(ig)/3.0);

		for(unsigned og=0; og<ngroups; og++){
			dir_ind[NDIMS+1] = og;
			global_index = sim->grid->scattering_delta[ID].direct_index(dir_ind);
			sim->grid->scattering_delta[ID][global_index] = tmp_delta[ig][og];
			sim->grid->scattering_phi0[ID][global_index] = tmp_phi0[ig][og] * pc::h;
		}
	}
}
