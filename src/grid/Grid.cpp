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
#include <iostream>
#include <iomanip>
#include <fstream>
#include "Grid.h"
#include "Lua.h"
#include "nulib_interface.h"
#include "global_options.h"
#include "Transport.h"
#include "H5Cpp.h"

using namespace std;
namespace pc = physical_constants;

//------------------------------------------------------------
// initialize the grid
//------------------------------------------------------------
void Grid::init(Lua* lua)
{
	// MPI stuff
	int my_rank,n_procs;
	MPI_Comm_rank( MPI_COMM_WORLD, &my_rank );
	MPI_Comm_size( MPI_COMM_WORLD, &n_procs);
	bool rank0 = (my_rank==0);

	// read the model file or fill in custom model
	read_model_file(lua);

	// read some parameters
	int do_relativity = lua->scalar<int>("do_relativity");
	int write_zones_every = lua->scalar<int>("write_zones_every");
	if(write_zones_every>0){
		output_zones_distribution = lua->scalar<int>("output_zones_distribution");
		output_hdf5 = lua->scalar<int>("output_hdf5");
	}
	do_annihilation = lua->scalar<int>("do_annihilation");

	// complain if the grid is obviously not right
	if(z.size()==0){
		cout << "Error: there are no grid zones." << endl;
		exit(5);
	}

	// calculate integrated quantities to check
	double total_rest_mass = 0.0;
	double total_KE        = 0.0;
	double total_TE        = 0.0;
	double total_hvis      = 0.0;
	double Tbar            = 0.0;
	double Yebar           = 0.0;
	int do_visc = lua->scalar<int>("do_visc");
    #pragma omp parallel for reduction(+:total_rest_mass,total_KE,total_TE,total_hvis,Tbar,Yebar)
	for(unsigned z_ind=0;z_ind<z.size();z_ind++){
		// zero out fluid velocity if not doing SR
		if(!do_relativity) for(unsigned i=0; i<ARRSIZE(z[z_ind].u); i++) z[z_ind].u[i] = 0;

		// calculate cell rest mass
		double rest_mass   = zone_rest_mass(z_ind);

		// sanity checks
		PRINT_ASSERT(z[z_ind].rho,>=,0.0);
		PRINT_ASSERT(zone_comoving_volume(z_ind),>=,0.0);
		PRINT_ASSERT(rest_mass,>=,0);

		// calculating totals and averages
		total_rest_mass += rest_mass;
		Tbar            += z[z_ind].T * rest_mass;
		Yebar           += z[z_ind].Ye * rest_mass;
		total_KE        += (Transport::lorentz_factor(z[z_ind].u,ARRSIZE(z[z_ind].u)) - 1.0) * rest_mass * pc::c*pc::c;
		total_TE        += rest_mass   / pc::m_n * pc::k * z[z_ind].T;
		if(do_visc) total_hvis += z[z_ind].H_vis * z[z_ind].rho * zone_lab_volume(z_ind);
	}

	// write out useful info about the grid
	if (rank0){
		cout << "#   mass = " << total_rest_mass/pc::M_sun << " Msun" << endl;
		cout << "#   <T> = " << Tbar/total_rest_mass*pc::k_MeV << " MeV" << endl;
		cout << "#   <Ye> = " << Yebar/total_rest_mass << endl;
		cout << "#   KE = " << total_KE << " erg" << endl;
		cout << "#   TE = " << total_TE << " erg" << endl;
		if(do_visc) cout << "#   hvis(nonrel) = " << total_hvis << " erg/s" << endl;
	}
}

//------------------------------------------------------------
// Write the grid information out to a file
//------------------------------------------------------------
void Grid::write_zones(const int iw) const
{
	PRINT_ASSERT(z.size(),>,0);

	// output all zone data in hdf5 format
	if(output_hdf5){
		PRINT_ASSERT(dimensionality(),>,0);
		string filename = Transport::filename("fluid",iw,".h5");
		H5::H5File file(filename, H5F_ACC_TRUNC);

		// write coordinates to the hdf5 file (implemented in each grid type)
		z[0].distribution[0].write_hdf5_coordinates(file);
		write_hdf5_coordinates(file);
		write_hdf5_data(file);
	}

	// output all zone data in text files
	else{
		ofstream outf;
		string filename = Transport::filename("fluid",iw,".dat");
		outf.open(filename.c_str());
		write_header(outf);
		for (unsigned z_ind=0; z_ind<z.size(); z_ind++){
			int dir_ind[dimensionality()];
			zone_directional_indices(z_ind, dir_ind, dimensionality());
			if(dimensionality()>0) if(dir_ind[dimensionality()-1]==0) outf << endl;
			write_line(outf,z_ind);
		}
		outf.close();
	}
}

double Grid::zone_rest_mass(const int z_ind) const{
	return z[z_ind].rho*zone_comoving_volume(z_ind);
}

double Grid::zone_comoving_volume(const int z_ind) const{
	// assumes v is orthonormal in cm/s
	if(!good_zone(z_ind)) return 0;
	else return zone_lab_volume(z_ind) * Transport::lorentz_factor(z[z_ind].u,3);
}

void Grid::write_header(ofstream& outf) const{
	outf << setprecision(4);
	outf << scientific;
	outf << "# ";
	double r[dimensionality()];
	zone_coordinates(0,r,dimensionality());
	unsigned c=0;
	for(unsigned i=0; i<dimensionality(); i++) outf << ++c << "-r[" << i << "]  ";
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
	if(do_annihilation) outf << ++c << "-annihilation_rate(erg/ccm/s,lab)  ";
	for(unsigned s=0; s<z[0].distribution.size(); s++) outf << ++c << "-e_rad"<< s << "(erg/ccm,lab) ";
	for(unsigned s=0; s<z[0].distribution.size(); s++) outf << ++c << "-avg_E"<< s << "(MeV,lab) ";
	if(output_zones_distribution){
		for(unsigned s=0; s<z[0].distribution.size(); s++){
			for(unsigned g=0; g<z[0].distribution[s].size(); g++){
				int nu_bin = z[0].distribution[s].nu_bin(g);
				int mu_bin = z[0].distribution[s].mu_bin(g);
				int phi_bin = z[0].distribution[s].phi_bin(g);
				outf << ++c << "-s"<<s<<"g"<<nu_bin<<"mu"<<mu_bin<<"phi"<<phi_bin<<"edens(erg/ccm)  ";
			}
		}
	}
	outf << ++c << "-Edens_invariant(erg/ccm/MeV^3) ";
	outf << endl;
}

void Grid::write_hdf5_data(H5::H5File file) const{
	// useful quantities
	H5::DataSet dataset;
	H5::DataSpace dataspace;
	vector<float> tmp(z.size(),0.0);

	// SET UP SCALAR DATASPACE
	hsize_t zdims[dimensionality()];
	dims(zdims,dimensionality());
	dataspace = H5::DataSpace(dimensionality(),&zdims[0]);

	// write comoving volume, assumes last index varies fastest.
	dataset = file.createDataSet("comoving_volume(ccm)",H5::PredType::IEEE_F32LE,dataspace);
	for(unsigned z_ind=0; z_ind<z.size(); z_ind++) tmp[z_ind] = zone_comoving_volume(z_ind);
	dataset.write(&tmp[0],H5::PredType::IEEE_F32LE);
	dataset.close();

	// write density
	dataset = file.createDataSet("rho(g|ccm,com)",H5::PredType::IEEE_F32LE,dataspace);
	for(unsigned z_ind=0; z_ind<z.size(); z_ind++) tmp[z_ind] = z[z_ind].rho;
	dataset.write(&tmp[0],H5::PredType::IEEE_F32LE);
	dataset.close();

	// write T_gas
	dataset = file.createDataSet("T_gas(MeV,com)",H5::PredType::IEEE_F32LE,dataspace);
	for(unsigned z_ind=0; z_ind<z.size(); z_ind++) tmp[z_ind] = z[z_ind].T*pc::k_MeV;
	dataset.write(&tmp[0],H5::PredType::IEEE_F32LE);
	dataset.close();

	// write Ye
	dataset = file.createDataSet("Ye",H5::PredType::IEEE_F32LE,dataspace);
	for(unsigned z_ind=0; z_ind<z.size(); z_ind++) tmp[z_ind] = z[z_ind].Ye;
	dataset.write(&tmp[0],H5::PredType::IEEE_F32LE);
	dataset.close();

	// write mue
	dataset = file.createDataSet("mue(MeV,com)",H5::PredType::IEEE_F32LE,dataspace);
	for(unsigned z_ind=0; z_ind<z.size(); z_ind++) tmp[z_ind] = nulib_eos_mue(z[z_ind].rho, z[z_ind].T, z[z_ind].Ye) * pc::ergs_to_MeV;
	dataset.write(&tmp[0],H5::PredType::IEEE_F32LE);
	dataset.close();

	// write H_vis
	dataset = file.createDataSet("H_vis(erg|s|g,com)",H5::PredType::IEEE_F32LE,dataspace);
	for(unsigned z_ind=0; z_ind<z.size(); z_ind++) tmp[z_ind] = z[z_ind].H_vis;
	dataset.write(&tmp[0],H5::PredType::IEEE_F32LE);
	dataset.close();

	// write H_abs
	dataset = file.createDataSet("H_abs(erg|s|g,com)",H5::PredType::IEEE_F32LE,dataspace);
	for(unsigned z_ind=0; z_ind<z.size(); z_ind++) tmp[z_ind] = z[z_ind].e_abs  / z[z_ind].rho;
	dataset.write(&tmp[0],H5::PredType::IEEE_F32LE);
	dataset.close();

	// write C_emit
	dataset = file.createDataSet("C_emit(erg|s|g,com)",H5::PredType::IEEE_F32LE,dataspace);
	for(unsigned z_ind=0; z_ind<z.size(); z_ind++) tmp[z_ind] = z[z_ind].e_emit / z[z_ind].rho;
	dataset.write(&tmp[0],H5::PredType::IEEE_F32LE);
	dataset.close();

	// write |v|
	dataset = file.createDataSet("|v|(cm|s,lab)",H5::PredType::IEEE_F32LE,dataspace);
	for(unsigned z_ind=0; z_ind<z.size(); z_ind++) tmp[z_ind] = sqrt(zone_speed2(z_ind));
	dataset.write(&tmp[0],H5::PredType::IEEE_F32LE);
	dataset.close();

	// write dYe_dt_abs
	dataset = file.createDataSet("dYe_dt_abs(1|s,lab)",H5::PredType::IEEE_F32LE,dataspace);
	for(unsigned z_ind=0; z_ind<z.size(); z_ind++){
		double n_baryons_per_ccm = z[z_ind].rho / Transport::mean_mass(z[z_ind].Ye);
		tmp[z_ind] = (z[z_ind].nue_abs-z[z_ind].anue_abs) / n_baryons_per_ccm / Transport::lorentz_factor(z[z_ind].u,ARRSIZE(z[z_ind].u));
	}
	dataset.write(&tmp[0],H5::PredType::IEEE_F32LE);
	dataset.close();

	// write annihilation_rate
	if(do_annihilation){
		dataset = file.createDataSet("annihilation_rate(erg|ccm|s,lab)",H5::PredType::IEEE_F32LE,dataspace);
		for(unsigned z_ind=0; z_ind<z.size(); z_ind++) tmp[z_ind] = z[z_ind].Q_annihil;
		dataset.write(&tmp[0],H5::PredType::IEEE_F32LE);
		dataset.close();
	}
	dataspace.close();

	// SET UP +1D DATASPACE
	hsize_t dims_plus1[dimensionality()+1];
	for(unsigned i=0; i<dimensionality(); i++) dims_plus1[i] = zdims[i];
	dims_plus1[dimensionality()] = z[0].distribution.size();
	dataspace = H5::DataSpace(dimensionality()+1,&dims_plus1[0]);
	tmp.resize(z.size() * z[0].distribution.size());

	// write neutrino energy density using contiguous temporary array
	dataset = file.createDataSet("e_rad(erg|ccm,lab)",H5::PredType::IEEE_F32LE,dataspace);
	for(unsigned z_ind=0; z_ind<z.size(); z_ind++)
		for(unsigned s=0; s<z[0].distribution.size(); s++){
			unsigned ind = z_ind*z[0].distribution.size() + s;
			tmp[ind] = z[z_ind].distribution[s].integrate();
		}
	dataset.write(&tmp[0],H5::PredType::IEEE_F32LE);
	dataset.close();

	// write average neutrino energy using contiguous temporary array
	dataset = file.createDataSet("avg_E(MeV,lab)",H5::PredType::IEEE_F32LE,dataspace);
	for(unsigned z_ind=0; z_ind<z.size(); z_ind++)
		for(unsigned s=0; s<z[0].distribution.size(); s++){
			unsigned ind = z_ind*z[0].distribution.size() + s;
			tmp[ind] = z[z_ind].distribution[s].average_nu() * pc::h_MeV;
		}
	dataset.write(&tmp[0],H5::PredType::IEEE_F32LE);
	dataset.close();
	dataspace.close();

	// write distribution function
	if(output_zones_distribution){

		// SET UP +4D DATASPACE
		vector<hsize_t> dims_plus4(dimensionality()+4,0);
		for(unsigned i=0; i<dimensionality(); i++) dims_plus4[i] = zdims[i];
		dims_plus4[dimensionality()  ] = z[0].distribution.size();
		dims_plus4[dimensionality()+1] = z[0].distribution[0].nu_dim();
		dims_plus4[dimensionality()+2] = z[0].distribution[0].mu_dim();
		dims_plus4[dimensionality()+3] =z[0].distribution[0].phi_dim();
		dataspace = H5::DataSpace(dimensionality()+4,&dims_plus4[0]);
		tmp.resize(z.size() * z[0].distribution.size() * z[0].distribution[0].nu_dim() * z[0].distribution[0].mu_dim() * z[0].distribution[0].phi_dim());
		dataset = file.createDataSet("distribution(erg|ccm,lab)",H5::PredType::IEEE_F32LE,dataspace);

		// set up subspace properties
		vector<hsize_t> offset(dims_plus4.size(),0);
		vector<hsize_t> stride(dims_plus4.size(),1);
		vector<hsize_t>  block(dims_plus4.size(),1);
		vector<hsize_t> spectrum_dims(dims_plus4);
		for(unsigned i=0; i<dimensionality()+1; i++) spectrum_dims[i] = 1;

		for(unsigned z_ind=0; z_ind<z.size(); z_ind++)
			for(unsigned s=0; s<z[0].distribution.size(); s++){
				// select the appropriate subspace
				int dir_ind[dimensionality()];
				zone_directional_indices(z_ind,dir_ind,dimensionality());
				for(unsigned i=0; i<dimensionality(); i++) offset[i] = dir_ind[i];
				offset[dimensionality()  ] = s;
				offset[dimensionality()+1] = 0;
				offset[dimensionality()+2] = 0;
				offset[dimensionality()+3] = 0;
				dataspace.selectHyperslab(H5S_SELECT_SET,&spectrum_dims[0],&offset[0],&stride[0],&block[0]);

				// write the data
				z[z_ind].distribution[s].write_hdf5_data(dataset, dataspace);
			}
		dataset.close();
	}
	dataspace.close();
}

void Grid::write_line(ofstream& outf, const int z_ind) const{
	double r[dimensionality()];
	zone_coordinates(z_ind,r,dimensionality());

	for(unsigned i=0; i<dimensionality(); i++) outf << r[i] << " ";

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

	double n_baryons_per_ccm = z[z_ind].rho / Transport::mean_mass(z[z_ind].Ye);
	double dYe_dt_abs = (z[z_ind].nue_abs-z[z_ind].anue_abs) / n_baryons_per_ccm / Transport::lorentz_factor(z[z_ind].u,ARRSIZE(z[z_ind].u));
	double dYe_dt_emit = -z[z_ind].l_emit / n_baryons_per_ccm / Transport::lorentz_factor(z[z_ind].u,ARRSIZE(z[z_ind].u));
	outf << dYe_dt_abs << "\t";
	outf << dYe_dt_emit << "\t";

	if(do_annihilation) outf << z[z_ind].Q_annihil << "\t";

	// integrated energy density and average energy for each species
	for(unsigned s=0; s<z[z_ind].distribution.size(); s++) outf << z[z_ind].distribution[s].integrate() << "\t";
	for(unsigned s=0; s<z[z_ind].distribution.size(); s++) outf << z[z_ind].distribution[s].average_nu() * pc::h_MeV << "\t";

	// integrated energy density
	//double integrated_edens = 0;
	//for(unsigned s=0; s<z[z_ind].distribution.size(); s++) integrated_edens += z[z_ind].distribution[s].integrate();
	//utf << integrated_edens << "\t";
	if(output_zones_distribution){
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

	double Binvar = 0;
	for(unsigned s=0; s<z[z_ind].distribution.size(); s++){
		for(unsigned i=0; i<z[z_ind].distribution[s].size(); i++){
			Binvar += z[z_ind].distribution[s].get(i) / pow(z[z_ind].distribution[s].nu_center(i) * pc::h_MeV,3);
		}
	}
	outf << Binvar << "\t";

	outf << endl;
}

bool Grid::good_zone(const int z_ind) const{
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
double Grid::zone_speed2(const int z_ind) const{
	PRINT_ASSERT(z_ind,>=,0);
	PRINT_ASSERT(z_ind,<,(int)z.size());
	double speed2 = Transport::dot(z[z_ind].u,z[z_ind].u,ARRSIZE(z[z_ind].u));
	return speed2;
}

double Grid::total_rest_mass() const{
	double mass = 0;
	#pragma omp parallel for reduction(+:mass)
	for(unsigned z_ind=0; z_ind<z.size(); z_ind++) mass += zone_rest_mass(z_ind);
	return mass;
}

