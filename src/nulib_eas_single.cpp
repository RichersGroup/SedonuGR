#include "global_options.h"
#include <iomanip>
#include <vector>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string.h>
#include <cmath>
#include "zone.h"
#include "nulib_interface.h"

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
	assert( nu_grid.size() > 0);
	assert( ye_grid.size() > 0);
	assert(  T_grid.size() > 0);
	assert(rho_grid.size() > 0);

	// read in the number of species and groups in the table
	cout << "# of species: "    << nulib_get_nspecies() << endl;
	cout << "# of groups: "     << nu_grid.size()       << " (" << nu_grid.min << "-" << nu_grid.x[nu_grid.size()-1] << ")" << endl;
	cout << "# of rho points: " << rho_grid.size()      << " (" << rho_grid[0] << "-" << rho_grid[rho_grid.size()-1] << ")" << endl;
	cout << "# of ye points: "  << ye_grid.size()       << " (" <<  ye_grid[0] << "-" <<  ye_grid[ ye_grid.size()-1] << ")" << endl;
	cout << "# of T points: "   << T_grid.size()        << " (" <<   T_grid[0] << "-" <<   T_grid[  T_grid.size()-1] << ")" << endl;

	//make vectors of appropriate sizes
	int n_groups = nu_grid.size();
	assert(n_groups > 0);
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
