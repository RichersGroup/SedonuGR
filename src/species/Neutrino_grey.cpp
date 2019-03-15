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
#include "nulib_interface.h"
#include "Transport.h"

using namespace std;
namespace pc = physical_constants;


Neutrino_grey::Neutrino_grey(){
	grey_opac = NaN;
	grey_abs_frac = NaN;
	grey_chempot = NaN;
}

//----------------------------------------------------------------
// called from species_general::init (neutrino-specific stuff)
//----------------------------------------------------------------
void Neutrino_grey::myInit(Lua* lua)
{
	grey_abs_frac = lua->scalar<double>("grey_abs_frac");
	grey_opac     = lua->scalar<double>("grey_opac");
	grey_chempot = lua->scalar<double>("grey_chempot") * pc::MeV_to_ergs;
}


//-----------------------------------------------------------------
// set emissivity, abs. opacity, and scat. opacity in zones
//-----------------------------------------------------------------
void Neutrino_grey::set_eas(const size_t z_ind, Grid* grid) const
{
	size_t ngroups = grid->nu_grid_axis.size();
	size_t dir_ind[NDIMS+1];
	grid->rho.indices(z_ind,dir_ind);

	PRINT_ASSERT(grey_abs_frac,>=,0);
	PRINT_ASSERT(grey_abs_frac,<=,1.0);
	grid->munue[z_ind] = grey_chempot;
	for(size_t j=0;j<grid->nu_grid_axis.size();j++)
	{
		dir_ind[NDIMS] = j;
		size_t global_index = grid->abs_opac[ID].direct_index(dir_ind);

		double a = grey_opac*grid->rho[z_ind]*grey_abs_frac;
		double s = grey_opac*grid->rho[z_ind]*(1.0-grey_abs_frac);
		PRINT_ASSERT(a,>=,0);
		PRINT_ASSERT(s,>=,0);

		grid->abs_opac[ID][global_index] = a;        // (1/cm)
		grid->scat_opac[ID][global_index] = s;        // (1/cm)

		grid->scattering_delta[ID].wipe();
		if(grid->scattering_delta[ID].size()>0){
		  for(size_t igout=0; igout<ngroups; igout++)
		    grid->partial_scat_opac[ID][igout][global_index] = (j==igout ? s : 0);
		}
	}
}

