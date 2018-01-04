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
#include "Neutrino_Nagakura.h"
#include "Transport.h"
#include "Grid.h"
#include "nulib_interface.h"
#include <fstream>
#include <sstream>
#include <string>

using namespace std;
namespace pc = physical_constants;


// constructor
Neutrino_Nagakura::Neutrino_Nagakura(){
}

//----------------------------------------------------------------
// called from species_general::init (neutrino-specific stuff)
//----------------------------------------------------------------
void Neutrino_Nagakura::myInit(Lua* lua)
{
    // set up the frequency table
    opacity_dir   = lua->scalar<string>("opacity_dir");

	// let the base class do the rest
	Neutrino::myInit(lua);
}


//-----------------------------------------------------------------
// set emissivity, abs. opacity, and scat. opacity in zones
//-----------------------------------------------------------------
void Neutrino_Nagakura::set_eas(int zone_index)
{
	unsigned dir_ind[NDIMS];
	sim->grid->rho.indices(zone_index,dir_ind);


	//=======================//
    // Open the opacity file //
	//=======================//
	stringstream filename;
	if(sim->grid->grid_type == "Grid1DSphere"){
		filename.str("");
	    filename << opacity_dir << "/opac_r" << zone_index << "_theta0.dat";
	}
	else if(sim->grid->grid_type == "Grid2DSphere"){
		vector<unsigned> dir_ind(2);
		hsize_t dims[2];
		sim->grid->dims(dims,2);
		sim->grid->zone_directional_indices(zone_index, dir_ind);
		filename.str("");
		filename << opacity_dir << "/opac_r" << dir_ind[0] << "_theta" << (dims[1]-dir_ind[1]-1) << ".dat"; // Hiroki's theta is backwards
	}
	else{
		cout << "ERROR: Incompatible grid and neutrino types. Hiroki only does spherical coordinates." << endl;
		cout << sim->grid->grid_type << endl;
		assert(false);
	}
    ifstream opac_file;
    opac_file.open(filename.str().c_str());

    //====================//
    // read the opacities //
    //====================//
    // ignore the first line
    string line;
    getline(opac_file,line);

    for(int inu=0; inu<nu_grid_axis->size(); inu++){
    	unsigned global_indices[NDIMS+1];
    	for(unsigned i=0; i<NDIMS; i++){
    		global_indices[i] = dir_ind[i];
    	}
    	global_indices[NDIMS] = inu;
    	unsigned global_index = sim->grid->abs_opac[ID].direct_index(global_indices);

            int itmp;
            double e=0, a=0, s=0;

            opac_file >> itmp; // group number
            PRINT_ASSERT(inu,==,itmp);

            double dtmp;
            if(ID==0){
                    // electron type
                    opac_file >> e; // emissivity (erg/ccm/s)
                    opac_file >> a;
                    opac_file >> s;

                    // anti-electron type
                    opac_file >> dtmp;
                    opac_file >> dtmp;
                    opac_file >> dtmp;

                    // heavy lepton type
                    opac_file >> dtmp;
                    opac_file >> dtmp;
                    opac_file >> dtmp;
            }
            else if(ID==1){
                    // electron type
                    opac_file >> dtmp;
                    opac_file >> dtmp;
                    opac_file >> dtmp;

                    // anti-electron type
                    opac_file >> e; // emissivity (erg/ccm/s)
                    opac_file >> a;
                    opac_file >> s;

                    // heavy lepton type
                    opac_file >> dtmp;
                    opac_file >> dtmp;
                    opac_file >> dtmp;
            }
            else if(ID==2){
                    // electron type
                    opac_file >> dtmp;
                    opac_file >> dtmp;
                    opac_file >> dtmp;

                    // electron type
                    opac_file >> dtmp;
                    opac_file >> dtmp;
                    opac_file >> dtmp;

                    // heavy lepton type
                    opac_file >> e; // emissivity (erg/ccm/s)
                    opac_file >> a;
                    opac_file >> s;
            }
            else{
                    cout << "ERROR: Neutrino ID not recognized!" << endl;
                    assert(false);
            }
            sim->grid->BB[ID][global_index] = e/(pc::h*nu_grid_axis->mid[inu]) * pc::c*pc::c/(4.*pc::pi * nu_grid_axis->delta3(inu)/3.0);
            sim->grid->abs_opac[ID][global_index] = a;
            sim->grid->scat_opac[ID][global_index] = s;
    }
    opac_file.close();
}
