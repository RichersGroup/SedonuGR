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
void Neutrino_grey::set_eas(int z_ind)
{
	unsigned ngroups = nu_grid_axis->size();

	PRINT_ASSERT(grey_abs_frac,>=,0);
	PRINT_ASSERT(grey_abs_frac,<=,1.0);
	for(unsigned j=0;j<nu_grid_axis->size();j++)
	{
		unsigned dir_ind[NDIMS+1];
		sim->grid->rho.indices(z_ind,dir_ind);
		dir_ind[NDIMS] = j;
		unsigned global_index = sim->grid->abs_opac[ID].direct_index(dir_ind);

		double nu  = sim->grid->nu_grid_axis.mid[j];        // (Hz)
		double bb  = Transport::number_blackbody(sim->grid->T[z_ind],0*pc::MeV_to_ergs,nu);

		double a = grey_opac*sim->grid->rho[z_ind]*grey_abs_frac;
		double s = grey_opac*sim->grid->rho[z_ind]*(1.0-grey_abs_frac);

		sim->grid->BB[s][global_index] = bb; // (#/s/cm^3/ster)
		sim->grid->abs_opac[ID][global_index] = a;        // (1/cm)
		sim->grid->scat_opac[ID][global_index] = s;        // (1/cm)
	}
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
