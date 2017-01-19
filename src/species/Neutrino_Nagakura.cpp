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
    string nugrid_filename = lua->scalar<string>("nugrid_filename");

    //==========================//
    // setup for frequency grid //
    //==========================//
    double minval = 0;
    double trash, tmp=0;
    vector<double> bintops = vector<double>(0);

    // get the frequency grid
    ifstream nugrid_file;
    nugrid_file.open(nugrid_filename.c_str());
    nugrid_file >> trash >> minval;
    minval /= pc::h_MeV;
    while(nugrid_file >> trash >> tmp){
            tmp /= pc::h_MeV;
            if(bintops.size()>0) PRINT_ASSERT(tmp,>,bintops[bintops.size()-1]);
            else PRINT_ASSERT(tmp,>,minval);
            bintops.push_back(tmp);
    }

    nugrid_file.close();

    // initialize the frequency grid
    nu_grid.init(bintops,minval,nu_grid.interpolation_method);

    //==================================//
    // set up the spectrum in each zone //
    //==================================//
	// temporary spectrum to be used for distribution function initialization
	LocateArray tmp_mugrid, tmp_phigrid;

	string mugrid_filename = lua->scalar<string>("distribution_mugrid_filename");
	ifstream mugrid_file;
	mugrid_file.open(mugrid_filename.c_str());
	mugrid_file >> trash >> minval;
	bintops = vector<double>(0);
	while(mugrid_file >> trash >> tmp){
		if(bintops.size()>0) PRINT_ASSERT(tmp,>,bintops[bintops.size()-1]);
		else PRINT_ASSERT(tmp,>,minval);
		bintops.push_back(tmp);
	}
	bintops[bintops.size()-1] = 1.0;
	PRINT_ASSERT(minval,==,-1);
	tmp_mugrid.init(bintops,minval);

	string phigrid_filename = lua->scalar<string>("distribution_phigrid_filename");
	ifstream phigrid_file;
	phigrid_file.open(phigrid_filename.c_str());
	phigrid_file >> trash >> minval;
	minval -= pc::pi;
	bintops = vector<double>(0);
	while(phigrid_file >> trash >> tmp){
		tmp -= pc::pi;
		if(bintops.size()>0) PRINT_ASSERT(tmp,>,bintops[bintops.size()-1]);
		else PRINT_ASSERT(tmp,>,minval);
		bintops.push_back(tmp);
	}
	bintops[bintops.size()-1] = pc::pi;
	PRINT_ASSERT(minval,==,-pc::pi);
	tmp_phigrid.init(bintops,minval);

	// use the above to initialize the zone's distribution function
	SpectrumArray tmp_spectrum;
	tmp_spectrum.init(nu_grid, tmp_mugrid, tmp_phigrid);

	//#pragma omp parallel for
	for(unsigned z_ind=0; z_ind<sim->grid->z.size(); z_ind++){
		sim->grid->z[z_ind].distribution.push_back(tmp_spectrum);
		PRINT_ASSERT(sim->grid->z[z_ind].distribution.size(),==,ID+1);
	}

	//===================================//
    // set neutrino's min and max values //
	//===================================//
    T_min   =  0.0;
    T_max   =  numeric_limits<double>::infinity();
    Ye_min  = 0;
    Ye_max  = 1;
    rho_min = -numeric_limits<double>::infinity();
    rho_max =  numeric_limits<double>::infinity();

	// let the base class do the rest
	Neutrino::myInit(lua);
}


//-----------------------------------------------------------------
// set emissivity, abs. opacity, and scat. opacity in zones
//-----------------------------------------------------------------
void Neutrino_Nagakura::set_eas(int zone_index)
{
	Zone* z = &(sim->grid->z[zone_index]);
	double ngroups = (double)emis[zone_index].size();

	//=======================//
    // Open the opacity file //
	//=======================//
	string filename;
	if(sim->grid->grid_type == "Grid1DSphere"){
	    filename = opacity_dir+string("/opac_r")+to_string(zone_index)+string(".dat");
	}
	else if(sim->grid->grid_type == "Grid2DSphere"){
		int dir_ind[2];
		sim->grid->zone_directional_indices(zone_index, dir_ind, 2);
		filename = opacity_dir+string("/opac_r")+to_string(dir_ind[0])+"_theta"+to_string(dir_ind[1])+string(".dat");
	}
	else{
		cout << "ERROR: Incompatible grid and neutrino types. Hiroki only does spherical coordinates." << endl;
		cout << sim->grid->grid_type << endl;
		assert(false);
	}
    ifstream opac_file;
    opac_file.open(filename);

    //====================//
    // read the opacities //
    //====================//
    // ignore the first line
    string line;
    getline(opac_file,line);

    for(int inu=0; inu<nu_grid.size(); inu++){
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
            emis[zone_index].set_value(inu,e);
            abs_opac[zone_index][inu] = a;
            scat_opac[zone_index][inu] = s;

            biased_emis[zone_index].set_value(inu,emis[zone_index].get_value(inu));
    }

    opac_file.close();


	emis[zone_index].normalize(cutoff/(double)ngroups);
	biased_emis[zone_index].normalize(cutoff/(double)ngroups);
}
