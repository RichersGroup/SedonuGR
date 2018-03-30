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
// do nothing
}


//-----------------------------------------------------------------
// set emissivity, abs. opacity, and scat. opacity in zones
//-----------------------------------------------------------------
void Neutrino_NuLib::set_eas(const unsigned z_ind, Grid* grid) const
{
	unsigned ngroups = grid->nu_grid_axis.size();
	unsigned dir_ind[NDIMS+2];
	grid->rho.indices(z_ind,dir_ind);

	vector<double> tmp_absopac(ngroups), tmp_scatopac(ngroups), tmp_BB(ngroups);
	vector< vector<double> > tmp_delta(ngroups, vector<double>(ngroups));
	vector<CDFArray> EoutCDF(ngroups);
	for(unsigned i=0; i<ngroups; i++) EoutCDF[i].resize(ngroups);
	nulib_get_eas_arrays(grid->rho[z_ind], grid->T[z_ind], grid->Ye[z_ind], ID,
			tmp_BB, tmp_absopac, tmp_scatopac, EoutCDF, tmp_delta);

	for(unsigned ig=0; ig<ngroups; ig++){
		dir_ind[NDIMS] = ig;
		unsigned global_index = grid->abs_opac[ID].direct_index(dir_ind);
		grid->abs_opac[ID][global_index] = tmp_absopac[ig];
		grid->scat_opac[ID][global_index] = tmp_scatopac[ig];
		grid->BB[ID][global_index] = tmp_BB[ig]; // erg/cm^2/s/sr - convert in next line
		grid->BB[ID][global_index] /= pc::h * pow(grid->nu_grid_axis.mid[ig],3) * grid->nu_grid_axis.delta(ig); // #/cm^2/s/sr/(Hz^3/3)

		if(grid->scattering_delta[ID].size()>0){
			PRINT_ASSERT(EoutCDF[ig].get(ngroups-1),==,1.0);
			PRINT_ASSERT(EoutCDF[ig].N,==,tmp_scatopac[ig]);
			for(unsigned og=0; og<ngroups; og++){
				dir_ind[NDIMS+1] = og;
				global_index = grid->scattering_delta[ID].direct_index(dir_ind);
				grid->scattering_delta[ID][global_index] = tmp_delta[ig][og];
				grid->scattering_EoutCDF[ID][global_index] = EoutCDF[ig].get(og);
			}
		}
	}
}

void Neutrino_NuLib::get_annihil_kernels(const double rho, const double T, const double Ye, const Axis& nuAxis, vector< vector< vector<double> > >& phi) const{
	nulib_get_epannihil_kernels(rho, T, Ye, ID, phi);
}
