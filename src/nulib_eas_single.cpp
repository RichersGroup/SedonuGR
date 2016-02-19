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
#include "nulib_interface.h"

namespace pc = physical_constants;

int main(int argc, char* argv[]){
	using namespace std;
	if(argc!=7){
		cout << "Usage: nulib_eas_single path_to_nulib_table.h5 rho(g/cm^3) T(MeV) Ye E(MeV) nulibID" << endl;
		exit(1);
	}

	//set test inputs
	double rho        = atof(argv[2]);           // g/cm^3
	double T          = atof(argv[3])/pc::k_MeV; // K
	double ye         = atof(argv[4]);
	double myenergy = atof(argv[5]);           // MeV
	double myfreq   = myenergy     /pc::h_MeV; // Hz
	int    nulibID  = atoi(argv[6]);

	//read in the nulib table
	cout << "initializing nulib" << endl;
	string filename = argv[1];
	nulib_init(filename);

	// grids
	locate_array nu_grid; // Hz
	vector<double> ye_grid;
	vector<double> T_grid; // K
	vector<double> rho_grid; // g/cm^3
	nulib_get_nu_grid(nu_grid);
	nulib_get_Ye_array(ye_grid);
	nulib_get_T_array(T_grid);
	nulib_get_rho_array(rho_grid);
	PRINT_ASSERT( nu_grid.size(),>,0);
	PRINT_ASSERT( ye_grid.size(),>,0);
	PRINT_ASSERT(  T_grid.size(),>,0);
	PRINT_ASSERT(rho_grid.size(),>,0);

	// read in the number of species and groups in the table
	cout << "# of species: "    << nulib_get_nspecies() << endl;
	cout << "# of groups: "     << nu_grid.size()       << " (" << nu_grid.min << "-" << nu_grid.x[nu_grid.size()-1] << ")" << endl;
	cout << "# of rho points: " << rho_grid.size()      << " (" << rho_grid[0] << "-" << rho_grid[rho_grid.size()-1] << ")" << endl;
	cout << "# of ye points: "  << ye_grid.size()       << " (" <<  ye_grid[0] << "-" <<  ye_grid[ ye_grid.size()-1] << ")" << endl;
	cout << "# of T points: "   << T_grid.size()        << " (" <<   T_grid[0] << "-" <<   T_grid[  T_grid.size()-1] << ")" << endl;

	//make vectors of appropriate sizes
	int n_groups = nu_grid.size();
	PRINT_ASSERT(n_groups,>,0);
	vector<double> absopac  (n_groups,0); // cm^-1
	vector<double> scatopac (n_groups,0); // cm^-1
	vector<double> pure_emis(n_groups,0); // erg/cm^3/s/ster/Hz
	cdf_array emis;                 // erg/cm^3/s/ster
	emis.resize(n_groups);

	//===================//
	// SINGLE LINE PLOTS //
	//===================//

	nulib_get_eas_arrays(rho, T, ye, nulibID, emis, absopac, scatopac);
	nulib_get_pure_emis (rho, T, ye, nulibID, pure_emis);
	cout << "e = " << nu_grid.value_at(myfreq, pure_emis) << " erg/cm^3/s/ster/Hz" << endl;
	cout << "a = " << nu_grid.value_at(myfreq, absopac)   << " 1/cm" << endl;
	cout << "s = " << nu_grid.value_at(myfreq, scatopac)  << " 1/cm" << endl;

	return 0;
}
