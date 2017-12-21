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
#include "Neutrino_GR1D.h"

using namespace std;
namespace pc = physical_constants;

//------------------------------------------------------------
// initialize the grid
//------------------------------------------------------------
void Grid::init(Lua* lua, Transport* insim)
{
	// MPI stuff
	int MPI_myID;
	MPI_Comm_rank( MPI_COMM_WORLD, &MPI_myID );
	bool rank0 = (MPI_myID==0);

	// set the transport pointer
	sim = insim;

	// read the model file or fill in custom model
	do_GR = lua->scalar<int>("do_GR");
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
		if(!do_relativity) for(unsigned i=0; i<3; i++) z[z_ind].u[i] = 0;

		// calculate cell rest mass
		double rest_mass   = zone_rest_mass(z_ind);

		// sanity checks
		PRINT_ASSERT(z[z_ind].rho,>=,0.0);
		PRINT_ASSERT(zone_com_3volume(z_ind),>=,0.0);
		PRINT_ASSERT(rest_mass,>=,0);

		// calculating totals and averages
		total_rest_mass += rest_mass;
		Tbar            += z[z_ind].T * rest_mass;
		Yebar           += z[z_ind].Ye * rest_mass;
		total_KE        += (LorentzHelper::lorentz_factor(z[z_ind].u,3) - 1.0) * rest_mass * pc::c*pc::c;
		total_TE        += rest_mass   / pc::m_n * pc::k * z[z_ind].T;
		if(do_visc) total_hvis += z[z_ind].H_vis * z[z_ind].rho * zone_com_3volume(z_ind);
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

    // get the frequency grid
	string neutrino_type = lua->scalar<string>("neutrino_type");
	if(neutrino_type=="NuLib") nulib_get_nu_grid(nu_grid_axis);
	else if(neutrino_type=="GR1D") Neutrino_GR1D::set_nu_grid(lua,&nu_grid_axis);
	else{
        double minval = 0;
        double trash, tmp=0;
        vector<double> bintops = vector<double>(0);
    	int n_nu   = lua->scalar<int>("nugrid_n");
    	if(n_nu>0){
    		double nu_start = lua->scalar<double>("nugrid_start");
    		double nu_stop  = lua->scalar<double>("nugrid_stop");
    		PRINT_ASSERT(nu_stop,>,nu_start);
    		nu_grid_axis = Axis(nu_start/pc::h_MeV,nu_stop/pc::h_MeV,n_nu);
    	}
    	else{
    		string nugrid_filename = lua->scalar<string>("nugrid_filename");
    		ifstream nugrid_file;
    		nugrid_file.open(nugrid_filename.c_str());
    		nugrid_file >> trash >> minval;
    		minval /= pc::h_MeV;
    		vector<double> binmid;
    		while(nugrid_file >> trash >> tmp){
    			tmp /= pc::h_MeV;
    			double last = bintops.size()>0 ? bintops[bintops.size()-1] : minval;
    			PRINT_ASSERT(tmp,>,last);
    			bintops.push_back(tmp);
    			binmid.push_back(0.5 * (last + tmp));
    		}
    		nugrid_file.close();
            nu_grid_axis = Axis(minval,bintops,binmid);
    	}
	}
	PRINT_ASSERT(nu_grid_axis.size(),>,0);
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
		z[0].distribution[0]->write_hdf5_coordinates(file, this);
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

// returning 0 causes the min distance to take over in propagate.cpp::which_event
double Grid::zone_cell_dist(const double x_up[3], const int z_ind) const{
	return 0;
}

double Grid::zone_rest_mass(const int z_ind) const{
	return z[z_ind].rho * zone_com_3volume(z_ind);
}

double Grid::zone_com_3volume(const int z_ind) const{
	// assumes v is orthonormal in cm/s
	if(z_ind<0) return 0;
	else return zone_lab_3volume(z_ind) * LorentzHelper::lorentz_factor(z[z_ind].u,3);
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
	outf << ++c << "-dYe_dt_abs(1/s,lab)  ";
	outf << ++c << "-dYe_dt_emit(1/s,lab)  ";
	if(do_annihilation) outf << ++c << "-annihilation_rate(erg/ccm/s,lab)  ";
	for(unsigned s=0; s<z[0].distribution.size(); s++) outf << ++c << "-e_rad"<< s << "(erg/ccm,lab) ";
	for(unsigned s=0; s<z[0].distribution.size(); s++) outf << ++c << "-avg_E"<< s << "(MeV,com) ";
	if(output_zones_distribution){
		for(unsigned s=0; s<z[0].distribution.size(); s++){
			outf << ++c << "-s"<<s;
			z[0].distribution[s]->write_header(outf);
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
	for(unsigned z_ind=0; z_ind<z.size(); z_ind++) tmp[z_ind] = zone_com_3volume(z_ind);
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

	// write dYe_dt_abs
	dataset = file.createDataSet("dYe_dt_abs(1|s,lab)",H5::PredType::IEEE_F32LE,dataspace);
	for(unsigned z_ind=0; z_ind<z.size(); z_ind++){
		double n_baryons_per_ccm = z[z_ind].rho / Transport::mean_mass(z[z_ind].Ye);
		tmp[z_ind] = (z[z_ind].nue_abs-z[z_ind].anue_abs) / n_baryons_per_ccm / LorentzHelper::lorentz_factor(z[z_ind].u,3);
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
			tmp[ind] = z[z_ind].distribution[s]->integrate();
		}
	dataset.write(&tmp[0],H5::PredType::IEEE_F32LE);
	dataset.close();

	// write average neutrino energy using contiguous temporary array
	dataset = file.createDataSet("avg_E(MeV,com)",H5::PredType::IEEE_F32LE,dataspace);
	for(unsigned z_ind=0; z_ind<z.size(); z_ind++)
		for(unsigned s=0; s<z[0].distribution.size(); s++){
			unsigned ind = z_ind*z[0].distribution.size() + s;
			tmp[ind] = z[z_ind].Edens_com[s] / z[z_ind].Ndens_com[s] * pc::h_MeV;
		}
	dataset.write(&tmp[0],H5::PredType::IEEE_F32LE);
	dataset.close();
	dataspace.close();

	// write distribution function
	if(output_zones_distribution){
		for(unsigned z_ind=0; z_ind<z.size(); z_ind++){
			int dir_ind[dimensionality()];
			zone_directional_indices(z_ind,dir_ind,dimensionality());
			for(unsigned s=0; s<z[0].distribution.size(); s++){
				z[z_ind].distribution[s]->write_hdf5_data(file, s, dir_ind, dimensionality());
			}
		}
	}
}

void Grid::write_line(ofstream& outf, const int z_ind) const{
	double r[dimensionality()];
	zone_coordinates(z_ind,r,dimensionality());

	for(unsigned i=0; i<dimensionality(); i++) outf << r[i] << " ";

	outf << zone_com_3volume(z_ind) << "\t";
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

	double n_baryons_per_ccm = z[z_ind].rho / Transport::mean_mass(z[z_ind].Ye);
	double dYe_dt_abs = (z[z_ind].nue_abs-z[z_ind].anue_abs) / n_baryons_per_ccm / LorentzHelper::lorentz_factor(z[z_ind].u,3);
	double dYe_dt_emit = -z[z_ind].l_emit / n_baryons_per_ccm / LorentzHelper::lorentz_factor(z[z_ind].u,3);
	outf << dYe_dt_abs << "\t";
	outf << dYe_dt_emit << "\t";

	if(do_annihilation) outf << z[z_ind].Q_annihil << "\t";

	// integrated energy density and average energy for each species
	for(unsigned s=0; s<z[z_ind].distribution.size(); s++) outf << z[z_ind].distribution[s]->integrate() << "\t";
	for(unsigned s=0; s<z[z_ind].distribution.size(); s++) outf << z[z_ind].Edens_com[s]/z[z_ind].Ndens_com[s] * pc::h_MeV << "\t";

	// integrated energy density
	//double integrated_edens = 0;
	//for(unsigned s=0; s<z[z_ind].distribution.size(); s++) integrated_edens += z[z_ind].distribution[s].integrate();
	//utf << integrated_edens << "\t";
	if(output_zones_distribution){
		for(unsigned s=0; s<z[0].distribution.size(); s++){
			z[z_ind].distribution[s]->write_line(outf);
		}
	}
//	for(unsigned s=0; s<z[z_ind].distribution.size(); s++){
//		vector<double> direction_integrated_edens;
//		z[z_ind].distribution[s].integrate_over_direction(direction_integrated_edens);
//		for(unsigned g=0; g<direction_integrated_edens.size(); g++){
//			outf << direction_integrated_edens[g] << "\t";
//		}
//	}

//	double Binvar = 0;
//	for(unsigned s=0; s<z[z_ind].distribution.size(); s++){
//		for(unsigned i=0; i<z[z_ind].distribution[s]->size(); i++){
//			Binvar += z[z_ind].distribution[s]->get(i) / pow(z[z_ind].distribution[s]->nu_center(i) * pc::h_MeV,3);
//		}
//	}
//	outf << Binvar << "\t";

	outf << endl;
}


double Grid::total_rest_mass() const{
	double mass = 0;
	#pragma omp parallel for reduction(+:mass)
	for(unsigned z_ind=0; z_ind<z.size(); z_ind++) mass += zone_rest_mass(z_ind);
	return mass;
}

// vector operations
template<int s>
double Grid::dot_Minkowski(const double a[], const double b[]){
	double product = 0;
	for(unsigned i=0; i<3; i++) product += a[i]*b[i];
	if(s==4) product -= a[3]*b[3];
	return product;
}
template double Grid::dot_Minkowski<3>(const double a[], const double b[]);
template double Grid::dot_Minkowski<4>(const double a[], const double b[]);

// normalize a vector
template<int s>
void Grid::normalize_Minkowski(double a[]){
	double inv_magnitude = 1./sqrt(abs( dot_Minkowski<s>(a,a) ));
	PRINT_ASSERT(inv_magnitude,<,INFINITY);
	for(unsigned i=0; i<s; i++) a[i] *= inv_magnitude;
}

void Grid::normalize_null_Minkowski(double a[4]){
	double spatial_norm = dot_Minkowski<3>(a,a);
	a[3] = sqrt(spatial_norm);
}

// radius given coordinates
double Grid::radius(const double x[4]){
	return sqrt(dot_Minkowski<3>(x,x));
}

void Grid::integrate_geodesic(LorentzHelper *lh) const{
	if(do_GR){
		const double* kold = lh->p_kup(lab);
		const double* xold = lh->p_xup();

		// get connection coefficients
		double gamma[4][4][4];
		connection_coefficients(xold,gamma);

		double lambda = lh->distance(lab) * (-n_dot(kold,xold));

		double dk_dlambda[4] = {0,0,0,0};
		double dx_dlambda[4] = {0,0,0,0};
		for(int a=0; a<4; a++){
			for(int mu=0; mu<4; mu++) for(int nu=0; nu<4; nu++)
				dk_dlambda[a] -= gamma[a][mu][nu] * kold[mu] * kold[nu];
			dx_dlambda[a]  = kold[a];
		}

		// get new x,k
		double knew[4],xnew[4];
		for(int i=0; i<4; i++){
			knew[i] = lh->p_kup(lab)[i] + dk_dlambda[i]*lambda;
			xnew[i] = lh->p_xup(   )[i] + dx_dlambda[i]*lambda;
		}

		// make vector null again
		normalize_null(knew,xnew);

		// actually change the values
		lh->set_p_xup(xnew,4);
		lh->set_p_kup<lab>(knew,4);
	}
	else{
		double xnew[4];
		for(int i=0; i<3; i++) xnew[i] = lh->p_xup()[i] + lh->distance(lab) * lh->p_kup(lab)[i] / lh->p_kup(lab)[3];
		xnew[3] = lh->p_xup()[3] + lh->distance(lab);
		lh->set_p_xup(xnew,4);
	}
}


// dot products
double Grid::dot(const double a[], const double b[], const double xup[], const int z_ind) const{
	double product = 0;
	if(do_GR){
		double g[4][4];
		g_down(xup,g,z_ind);
		for(int mu=0; mu<4; mu++) for(int nu=0; nu<4; nu++)
			product += a[mu] * b[nu] * g[mu][nu];
	}
	else product = dot_Minkowski<4>(a,b);
	return product;
}

double Grid::dot3(const double a[], const double b[], const double xup[], const int z_ind) const{
	double product = 0;
	if(do_GR){
		double g[4][4];
		g_down(xup,g,z_ind);
		for(int mu=0; mu<3; mu++) for(int nu=0; nu<3; nu++)
			product += a[mu] * b[nu] * g[mu][nu];
	}
	else product = dot_Minkowski<3>(a,b);
	return product;
}
void Grid::normalize(double a[], const double xup[], const int z_ind) const{
	double inv_norm = 0;
	if(do_GR){
		double g[4][4];
		g_down(xup,g,z_ind);
		for(int mu=0; mu<4; mu++) for(int nu=0; nu<4; nu++)
			inv_norm += a[mu] * a[nu] * g[mu][nu];
		inv_norm = 1.0 / sqrt(inv_norm);

		for(int mu=0; mu<4; mu++) a[mu] *= inv_norm;
	}
	else normalize_Minkowski<4>(a);
}
void Grid::normalize_null(double a[], const double xup[], const int z_ind) const{
	if(do_GR){
		double g[4][4];
		g_down(xup,g,z_ind);
		double A = g[3][3];
		double B=0, C=0;

		for(int i=0; i<3; i++){
			B += a[i] * g[i][3] / A;

			for(int j=0; j<3; j++)
				C += a[i] * a[j] * g[i][j] / A;
		}

		double result = 0.5 * (B + sqrt(B*B - 4.*C));
		PRINT_ASSERT(result,>,0);

		a[3] = result;
	}
	else normalize_null_Minkowski(a);
}

// isotropic scatter, done in COMOVING frame
void Grid::isotropic_direction(double D[3], ThreadRNG *rangen) const
{
	// Randomly generate new direction isotropically in comoving frame
	double mu  = 1 - 2.0*rangen->uniform();
	double phi = 2.0*pc::pi*rangen->uniform();
	double smu = sqrt(1 - mu*mu);

	D[0] = smu*cos(phi);
	D[1] = smu*sin(phi);
	D[2] = mu;
	Grid::normalize_Minkowski<3>(D);
}

void Grid::isotropic_kup_tet(const double nu, double kup_tet[4], const double xup[4], ThreadRNG *rangen) const{
	double D[3];
	isotropic_direction(D,rangen);

	kup_tet[0] = nu * D[0] * 2.0*pc::pi / pc::c;
	kup_tet[1] = nu * D[1] * 2.0*pc::pi / pc::c;
	kup_tet[2] = nu * D[2] * 2.0*pc::pi / pc::c;
	kup_tet[3] = nu        * 2.0*pc::pi / pc::c;

	if(do_GR){
		// move from tetrad frame to comoving frame
		double v[3];
		interpolate_fluid_velocity(xup,v);
		double u[4];
		double gamma = LorentzHelper::lorentz_factor(v,3);
		u[0] = gamma * v[0];
		u[1] = gamma * v[1];
		u[2] = gamma * v[2];
		u[3] = gamma * pc::c;
	}
}

void Grid::orthogonalize(double v[4], const double e[4], const double xup[4]) const{
	PRINT_ASSERT(abs(dot(e,e,xup)-1.0),<,TINY); // assume the basis vector is normalized
	PRINT_ASSERT(do_GR,==,true);

	double projection = dot(v,e,xup);
	for(int mu=0; mu<4; mu++) v[mu] -= projection * v[mu];
}
void Grid::tetrad_basis(const double xup[4], const double u[4], double e[4][4]) const{
	// normalize four-velocity to get timelike vector
	for(int mu=0; mu<4; mu++) e[3][mu] = u[mu];
	normalize(e[3], xup);

	// use x0 as a trial vector
	e[0][0] = 1.0;
	e[0][1] = 0;
	e[0][2] = 0;
	e[0][3] = 0;
	orthogonalize(e[0],e[3],xup);
	normalize(e[0],xup);

	// use x1 as a trial vector
	e[1][0] = 0;
	e[1][1] = 1.0;
	e[1][2] = 0;
	e[1][3] = 0;
	orthogonalize(e[1],e[3],xup);
	orthogonalize(e[1],e[0],xup);
	normalize(e[1],xup);

	// use x2 as a trial vector
	e[2][0] = 0;
	e[2][1] = 0;
	e[2][2] = 1.0;
	e[2][3] = 0;
	orthogonalize(e[2],e[3],xup);
	orthogonalize(e[2],e[0],xup);
	orthogonalize(e[2],e[1],xup);
	normalize(e[2],xup);

	// sanity checks
	PRINT_ASSERT(abs(dot(e[0],e[1],xup)),<,TINY);
	PRINT_ASSERT(abs(dot(e[0],e[2],xup)),<,TINY);
	PRINT_ASSERT(abs(dot(e[0],e[3],xup)),<,TINY);
	PRINT_ASSERT(abs(dot(e[1],e[2],xup)),<,TINY);
	PRINT_ASSERT(abs(dot(e[1],e[3],xup)),<,TINY);
	PRINT_ASSERT(abs(dot(e[2],e[3],xup)),<,TINY);
}
void Grid::coord_to_tetrad(const double xup[4], const double u[4], const double kup_coord[4], double kup_tet[4]) const{
	double e[4][4];
	tetrad_basis(xup, u, e);

	// transform to coordinate frame
	for(int mu=0; mu<4; mu++) kup_tet[mu] = dot(kup_coord,e[mu],xup);
}
void Grid::tetrad_to_coord(const double xup[4], const double u[4], const double kup_tet[4], double kup_coord[4]) const{
	double e[4][4];
	tetrad_basis(xup, u, e);

	// transform to coordinate frame
	for(int mu=0; mu<4; mu++){
		kup_coord[mu] = 0;
		for(int nu=0; nu<4; nu++)
			kup_coord[mu] += kup_tet[nu] * e[nu][mu];
	}
}

void Grid::random_core_x_D(const double r_core, ThreadRNG *rangen, double x3[3], double D[3]) const{
	PRINT_ASSERT(do_GR,==,false);

	double phi_core   = 2*pc::pi*rangen->uniform();
	double cosp_core  = cos(phi_core);
	double sinp_core  = sin(phi_core);
	double cost_core  = 1 - 2.0*rangen->uniform();
	double sint_core  = sqrt(1-cost_core*cost_core);

	double a_phot = r_core + r_core*1e-10;
	x3[0] = a_phot*sint_core*cosp_core;
	x3[1] = a_phot*sint_core*sinp_core;
	x3[2] = a_phot*cost_core;

	// pick photon propagation direction wtr to local normal
	double phi_loc = 2*pc::pi*rangen->uniform();
	// choose sqrt(R) to get outward, cos(theta) emission
	double cost_loc  = sqrt(rangen->uniform());
	double sint_loc  = sqrt(1 - cost_loc*cost_loc);
	// local direction vector
	double D_xl = sint_loc*cos(phi_loc);
	double D_yl = sint_loc*sin(phi_loc);
	double D_zl = cost_loc;
	// apply rotation matrix to convert D vector into overall frame
	D[0] = cost_core*cosp_core*D_xl-sinp_core*D_yl+sint_core*cosp_core*D_zl;
	D[1] = cost_core*sinp_core*D_xl+cosp_core*D_yl+sint_core*sinp_core*D_zl;
	D[2] = -sint_core*D_xl+cost_core*D_zl;
	Grid::normalize_Minkowski<3>(D);
}

double Grid::n_dot(const double invec[4], const double xup[4], const int z_ind) const{
	return -lapse(xup,z_ind) * invec[3];
}

void Grid::g_down(const double xup[4], double g[4][4], const int z_ind) const{
	g3_down(xup,g,z_ind);

	double alpha = lapse(xup,z_ind);

	double betaup[4];
	shiftup(betaup,xup,z_ind);

	// betaup[3] is zero by definition
	double beta_dot_beta=0;
	for(int i=0; i<3; i++){
		double betadowni = 0;
		for(int j=0; j<3; j++) betadowni += betaup[j] * g[i][j];
		beta_dot_beta += betaup[i] * betadowni;
		g[3][i] = g[i][3] = betadowni;
	}
	g[3][3] = -alpha*alpha + beta_dot_beta;
	PRINT_ASSERT(g[3][3],<,0);
}

double Grid::zone_4volume(const int z_ind) const{
	return zone_lab_3volume(z_ind) * zone_lapse(z_ind);
}
