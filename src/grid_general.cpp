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

#include <mpi.h>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
#include "grid_general.h"
#include "Lua.h"
#include "nulib_interface.h"
#include "global_options.h"
#include "H5Cpp.h"

//------------------------------------------------------------
// initialize the grid
//------------------------------------------------------------
void grid_general::init(Lua* lua)
{
	// MPI stuff
	int my_rank,n_procs;
	MPI_Comm_rank( MPI_COMM_WORLD, &my_rank );
	MPI_Comm_size( MPI_COMM_WORLD, &n_procs);
	bool rank0 = (my_rank==0);

	// read the model file or fill in custom model
	std::string model_file = lua->scalar<std::string>("model_file");
	int do_relativity = lua->scalar<int>("do_relativity");
	if(model_file == "custom") custom_model(lua);
	else read_model_file(lua);
	output_distribution = lua->scalar<int>("output_distribution");
	output_hdf5 = lua->scalar<int>("output_hdf5");

	// complain if the grid is obviously not right
	if(z.size()==0){
		cout << "Error: there are no grid zones." << endl;
		exit(5);
	}

	// calculate integrated quantities to check
	double total_nonrel_mass = 0.0;
	double total_rest_mass   = 0.0;
	double total_rel_KE      = 0.0;
	double total_nonrel_KE   = 0.0;
	double total_rel_TE      = 0.0;
	double total_nonrel_TE   = 0.0;
	double total_hvis        = 0.0;
	double nonrel_Tbar       = 0.0;
	double rel_Tbar          = 0.0;
	double nonrel_Yebar      = 0.0;
	double rel_Yebar         = 0.0;
	int do_visc = lua->scalar<int>("do_visc");
    #pragma omp parallel for reduction(+:total_nonrel_mass,total_rest_mass,total_rel_KE,total_nonrel_KE,total_rel_TE,total_nonrel_TE,total_hvis,nonrel_Tbar,rel_Tbar,nonrel_Yebar,rel_Yebar)
	for(unsigned z_ind=0;z_ind<z.size();z_ind++){
		double rest_mass   = zone_rest_mass(z_ind);
		assert(rest_mass >= 0);
		double nonrel_mass = z[z_ind].rho * zone_lab_volume(z_ind);
		assert(nonrel_mass >= 0);
		vector<double> r;
		zone_coordinates(z_ind,r);

		//if(grid->z[z_ind].rho > 1.0e8){ // && r[1] > pc::pi/3.0 && r[1] < pc::pi/2.0){
		total_rest_mass += rest_mass;
		total_nonrel_mass += nonrel_mass;
		nonrel_Tbar += z[z_ind].T * rest_mass;
		rel_Tbar += z[z_ind].T * nonrel_mass;
		nonrel_Yebar += z[z_ind].Ye * rest_mass;
		rel_Yebar += z[z_ind].Ye * nonrel_mass;
		total_rel_KE    += (rest_mass>0 ? (transport::lorentz_factor(z[z_ind].v) - 1.0) * rest_mass * pc::c*pc::c : 0);
		total_nonrel_KE += 0.5 * nonrel_mass * zone_speed2(z_ind);
		total_rel_TE    += (rest_mass>0 ? rest_mass   / pc::m_n * pc::k * z[z_ind].T : 0);
		total_nonrel_TE += nonrel_mass / pc::m_n * pc::k * z[z_ind].T;
		if(do_visc) total_hvis += z[z_ind].H_vis * z[z_ind].rho * zone_comoving_volume(z_ind);
		if(!do_relativity) for(unsigned i=0; i<z[z_ind].v.size(); i++) z[z_ind].v[i] = 0;
		//}
	}
	if (rank0){
		cout << "#   mass = " << total_rest_mass/pc::M_sun << " M_sun (nonrel: " << total_nonrel_mass/pc::M_sun << " M_sun)" << endl;
		cout << "#   <T> = " << rel_Tbar/total_rest_mass*pc::k_MeV << " MeV (nonrel: " << nonrel_Tbar/total_nonrel_mass*pc::k_MeV << " MeV)" <<endl;
		cout << "#   <Ye> = " << rel_Yebar/total_rest_mass << " (nonrel: " << nonrel_Yebar/total_nonrel_mass << ")" <<endl;
		cout << "#   KE = " << total_rel_KE << " erg (nonrel: " << total_nonrel_KE << " erg)" << endl;
		cout << "#   TE = " << total_rel_TE << " erg (nonrel: " << total_nonrel_TE << " erg)" << endl;
		if(do_visc) cout << "#   hvis = " << total_hvis << " erg/s" << endl;
	}
}


//------------------------------------------------------------
// Write the grid information out to a file
//------------------------------------------------------------
void grid_general::write_zones(const int iw) const
{
	assert(z.size()>0);

	// output all zone data in hdf5 format
	if(output_hdf5){
		H5::H5File file();
	}

	// output all zone data in text files
	else{
		ofstream outf;
		string filename = transport::filename("fluid",iw,outf);
		outf.open(filename.c_str());
		write_header(outf);
		vector<int> dir_ind;
		for (unsigned z_ind=0; z_ind<z.size(); z_ind++)
		{
			zone_directional_indices(z_ind, dir_ind);
			if(dir_ind.size()>0) if(dir_ind[dir_ind.size()-1]==0) outf << endl;
			write_line(outf,z_ind);
		}
		outf.close();
	}
}

double grid_general::zone_rest_mass(const int z_ind) const{
	return z[z_ind].rho*zone_comoving_volume(z_ind);
}

double grid_general::zone_comoving_volume(const int z_ind) const{
	// assumes v is orthonormal in cm/s
	vector<double> v_tmp(3,0);
	for(unsigned i=0; i<z[z_ind].v.size(); i++) v_tmp[i] = z[z_ind].v[i];
	if(!good_zone(z_ind)) return 0;
	else return zone_lab_volume(z_ind) * transport::lorentz_factor(v_tmp);
}

void grid_general::write_header(ofstream& outf) const{
	outf << setprecision(4);
	outf << scientific;
	outf << "# ";
	vector<double> r;
	zone_coordinates(0,r);
	unsigned c=0;
	for(unsigned i=0; i<r.size(); i++) outf << ++c << "-r[" << i << "]  ";
	outf << ++c << "-comoving_volume(ccm)  ";
	outf << ++c << "-rho(g/ccm,com)  ";
	outf << ++c << "-T_gas(MeV,com)  ";
	outf << ++c << "-Ye  ";
	outf << ++c << "-mue(MeV,com)  ";
	outf << ++c << "-H_vis(erg/s/g,com)  ";
	outf << ++c << "-H_abs(erg/g/s,com)  ";
	outf << ++c << "-C_emit(erg/g/s,com)  ";
	outf << ++c << "-|v|(cm/s,lab)  ";
	outf << ++c << "-dYe_dt_abs(1/s,lab)  ";
	outf << ++c << "-dYe_dt_emit(1/s,lab)  ";
	outf << ++c << "-annihilation_rate(erg/ccm/s,lab)  ";
	for(unsigned s=0; s<z[0].distribution.size(); s++) outf << ++c << "-e_rad"<< s << "(erg/ccm,lab) ";
	for(unsigned s=0; s<z[0].distribution.size(); s++) outf << ++c << "-avg_E"<< s << "(MeV,lab) ";
	if(output_distribution){
		for(unsigned s=0; s<z[0].distribution.size(); s++){
			for(unsigned g=0; g<z[0].distribution[s].size(); g++){
				int nu_bin = z[0].distribution[s].nu_bin(g);
				int mu_bin = z[0].distribution[s].mu_bin(g);
				int phi_bin = z[0].distribution[s].phi_bin(g);
				outf << ++c << "-s"<<s<<"g"<<nu_bin<<"mu"<<mu_bin<<"phi"<<phi_bin<<"edens(erg/ccm)  ";
			}
		}
	}
	outf << endl;
}

void grid_general::write_line(ofstream& outf, const int z_ind) const{
	vector<double> r;
	zone_coordinates(z_ind,r);

	for(unsigned i=0; i<r.size(); i++) outf << r[i] << " ";

	outf << zone_comoving_volume(z_ind) << "\t";
	outf << z[z_ind].rho   << "\t";
	outf << z[z_ind].T*pc::k_MeV << "\t";
	outf << z[z_ind].Ye    << "\t";
	if(nulib_in_range(z[z_ind].rho, z[z_ind].T, z[z_ind].Ye))
		outf << nulib_eos_mue(z[z_ind].rho, z[z_ind].T, z[z_ind].Ye) * pc::ergs_to_MeV << "\t";
	else outf << "UNDEF\t";
	outf << z[z_ind].H_vis << "\t";

	double H_abs  = z[z_ind].e_abs  / z[z_ind].rho;
	double H_emit = z[z_ind].e_emit / z[z_ind].rho;
	outf << H_abs << "\t";
	outf << H_emit << "\t";

	outf << sqrt(zone_speed2(z_ind)) << "\t";

	double n_baryons_per_ccm = z[z_ind].rho / transport::mean_mass(z[z_ind].Ye);
	double dYe_dt_abs = (z[z_ind].nue_abs-z[z_ind].anue_abs) / n_baryons_per_ccm / transport::lorentz_factor(z[z_ind].v);
	double dYe_dt_emit = -z[z_ind].l_emit / n_baryons_per_ccm / transport::lorentz_factor(z[z_ind].v);
	outf << dYe_dt_abs << "\t";
	outf << dYe_dt_emit << "\t";

	outf << z[z_ind].Q_annihil << "\t";

	// integrated energy density and average energy for each species
	for(unsigned s=0; s<z[z_ind].distribution.size(); s++) outf << z[z_ind].distribution[s].integrate() << "\t";
	for(unsigned s=0; s<z[z_ind].distribution.size(); s++) outf << z[z_ind].distribution[s].average_nu() * pc::h_MeV << "\t";

	// integrated energy density
	//double integrated_edens = 0;
	//for(unsigned s=0; s<z[z_ind].distribution.size(); s++) integrated_edens += z[z_ind].distribution[s].integrate();
	//utf << integrated_edens << "\t";
	if(output_distribution){
		for(unsigned s=0; s<z[0].distribution.size(); s++){
			for(unsigned g=0; g<z[0].distribution[s].size(); g++){
				outf << z[z_ind].distribution[s].get(g) << "\t";
			}
		}
	}
//	for(unsigned s=0; s<z[z_ind].distribution.size(); s++){
//		vector<double> direction_integrated_edens;
//		z[z_ind].distribution[s].integrate_over_direction(direction_integrated_edens);
//		for(unsigned g=0; g<direction_integrated_edens.size(); g++){
//			outf << direction_integrated_edens[g] << "\t";
//		}
//	}

	outf << endl;
}

bool grid_general::good_zone(const int z_ind) const{
	//const zone* z = &(grid->z[z_ind]);
	//return (z->rho >= rho_min && z->rho <= rho_max &&
	//  	    z->Ye  >= Ye_min  && z->Ye  <=  Ye_max &&
	//        z->T   >=  T_min  && z->T   <=   T_max);
	if(z_ind < 0) return false;
	else if(zone_speed2(z_ind) > pc::c*pc::c) return false;
	else return true;
}



//------------------------------------
// get the velocity squared of a zone
//------------------------------------
double grid_general::zone_speed2(const int z_ind) const{
	assert(z_ind >= 0);
	assert(z_ind < (int)z.size());
	double speed2 = transport::dot(z[z_ind].v,z[z_ind].v);
	return speed2;
}

double grid_general::total_rest_mass() const{
	double mass = 0;
	#pragma omp parallel for reduction(+:mass)
	for(unsigned z_ind=0; z_ind<z.size(); z_ind++) mass += zone_rest_mass(z_ind);
	return mass;
}
