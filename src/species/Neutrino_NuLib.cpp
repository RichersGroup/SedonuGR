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
void Neutrino_NuLib::myInit(Lua* /*lua*/)
{
// do nothing
}


//-----------------------------------------------------------------
// set emissivity, abs. opacity, and scat. opacity in zones
//-----------------------------------------------------------------
void Neutrino_NuLib::set_eas(const size_t z_ind, Grid* grid) const
{
	size_t ngroups = grid->nu_grid_axis.size();
	size_t dir_ind[NDIMS+2];
	grid->rho.indices(z_ind,dir_ind);

	vector<double> tmp_absopac(ngroups), tmp_scatopac(ngroups);
	vector< vector<double> > tmp_delta(ngroups, vector<double>(ngroups)); //[igin][igout]
	vector< vector<double> > tmp_partial_opac(ngroups, vector<double>(ngroups));  //[igin][igout]
	nulib_get_eas_arrays(grid->rho[z_ind], grid->T[z_ind], grid->Ye[z_ind], ID,
			tmp_absopac, tmp_scatopac, tmp_partial_opac, tmp_delta);

	if(ID==0) grid->munue[z_ind] = nulib_eos_munue(grid->rho[z_ind], grid->T[z_ind], grid->Ye[z_ind]);
	for(size_t igin=0; igin<ngroups; igin++){
		dir_ind[NDIMS] = igin;
		size_t global_index1 = grid->abs_opac[ID].direct_index(dir_ind);
		grid->abs_opac[ID][global_index1] = tmp_absopac[igin];
		grid->scat_opac[ID][global_index1] = tmp_scatopac[igin];

		if(grid->scattering_delta[ID].size()>0){
			for(size_t igout=0; igout<ngroups; igout++){
				grid->partial_scat_opac[ID][igout][global_index1] = tmp_partial_opac[igin][igout];
				dir_ind[NDIMS+1] = igout;
				size_t global_index2 = grid->scattering_delta[ID].direct_index(dir_ind);
				grid->scattering_delta[ID][global_index2] = tmp_delta[igin][igout];
			}
		}
	}
}

void Neutrino_NuLib::get_annihil_kernels(const double rho, const double T, const double Ye, const Axis& /*nuAxis*/, vector< vector< vector<double> > >& phi) const{
	nulib_get_epannihil_kernels(rho, T, Ye, ID, phi);
}
