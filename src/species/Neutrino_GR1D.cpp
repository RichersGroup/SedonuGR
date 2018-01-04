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
#include "Neutrino_GR1D.h"
#include "Transport.h"
#include "Grid.h"
#include "nulib_interface.h"
#include "H5Cpp.h"

using namespace std;
namespace pc = physical_constants;


// constructor
Neutrino_GR1D::Neutrino_GR1D(){
	ghosts1 = -1;
	n_GR1D_zones = -1;
	GR1D_emit_outside_shock = 1;
}

//----------------------------------------------------------------
// called from species_general::init (neutrino-specific stuff)
//----------------------------------------------------------------
void Neutrino_GR1D::myInit(Lua* lua)
{
	// emit outside shock?
	GR1D_emit_outside_shock = lua->scalar<int>("GR1D_emit_outside_shock");

	// let the base class do the rest
	Neutrino::myInit(lua);
}

void Neutrino_GR1D::set_nu_grid(Lua* lua, Axis* nu_grid){
	// get nugrid from nulib file
	string nulibfilename = lua->scalar<string>("nulib_file");
	H5::H5File file( nulibfilename.c_str(), H5F_ACC_RDONLY );

	// read tops from hdf5 file
	H5::DataSet dataset = file.openDataSet("bin_top");
	H5::DataSpace dataspace = dataset.getSpace();
	PRINT_ASSERT(dataspace.getSimpleExtentNdims(),==,1);
	hsize_t dims[1];
	dataspace.getSimpleExtentDims(dims);
	vector<double> tops(dims[0]);
	dataset.read(&(tops[0]),H5::PredType::IEEE_F64LE);

	// read the centers from hdf5 file
	dataset = file.openDataSet("energies");
	dataspace = dataset.getSpace();
	PRINT_ASSERT(dataspace.getSimpleExtentNdims(),==,1);
	dataspace.getSimpleExtentDims(dims);
	vector<double> mid(dims[0]);
	dataset.read(&(mid[0]),H5::PredType::IEEE_F64LE);

	dataset.close();
	file.close();

	// adjust units
	for(int i=0; i<nu_grid->size(); i++){
		tops[i] /= pc::h_MeV;
		mid[i] /= pc::h_MeV;
	}
	*nu_grid = Axis(0, tops, mid);
}

//-----------------------------------------------------------------
// set emissivity, abs. opacity, and scat. opacity in zones
//-----------------------------------------------------------------
void Neutrino_GR1D::set_eas(int zone_index)
{
	// do nothing - opacities are set externally
}

int extract_MC_index(const int z_ind, const int ID, const int nspecies, const int inu, const int ngroups){
	return inu + ID*ngroups + z_ind*ngroups*nspecies;
}
void Neutrino_GR1D::set_eas_external(const double* easarray, const double GR1D_tau_crit, bool* extract_MC, const double rshock){
	PRINT_ASSERT(ghosts1,>=,0);
	PRINT_ASSERT(n_GR1D_zones,>=,sim->grid->rho.size());
	const int ngroups=nu_grid_axis->size();
	const int nspecies = sim->species_list.size();

	// first, get opacities everywhere
	for(int z_ind=0; z_ind<sim->grid->rho.size(); z_ind++){
		for(int inu=0; inu<ngroups; inu++){
			unsigned dir_ind[NDIMS+1];
			sim->grid->rho.indices(z_ind,dir_ind);
			dir_ind[NDIMS] = inu;
			unsigned global_index = sim->grid->abs_opac[ID].direct_index(dir_ind);

			// indexed as eas(zone,species,group,e/a/s). The leftmost one varies fastest.
			int aind = (z_ind+ghosts1) + ID*n_GR1D_zones + inu*nspecies*n_GR1D_zones + 1*ngroups*nspecies*n_GR1D_zones;
			int sind = (z_ind+ghosts1) + ID*n_GR1D_zones + inu*nspecies*n_GR1D_zones + 2*ngroups*nspecies*n_GR1D_zones;
			int eind = (z_ind+ghosts1) + ID*n_GR1D_zones + inu*nspecies*n_GR1D_zones + 0*ngroups*nspecies*n_GR1D_zones;

			PRINT_ASSERT(easarray[aind],>=,0);
			PRINT_ASSERT(easarray[sind],>=,0);

			// set opacities
			sim->grid->abs_opac[ID][global_index] = easarray[aind] / nulib_opacity_gf;
			sim->grid->abs_opac[ID][global_index] = easarray[sind] / nulib_opacity_gf;
			sim->grid->BB[ID][global_index] = easarray[eind] / nulib_emissivity_gf /(pc::h*nu_grid_axis->mid[inu]) * pc::c*pc::c/(4.*pc::pi * nu_grid_axis->delta3(inu)/3.0);
		}
	}
}


