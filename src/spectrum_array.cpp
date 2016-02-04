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

#include <limits>
#include <omp.h>
#include <mpi.h>
#include <iostream>
#include <sstream>
#include <math.h>
#include <vector>
#include <fstream>
#include <string>
#include "global_options.h"
#include "spectrum_array.h"

//--------------------------------------------------------------
// Initialization and Allocation
//--------------------------------------------------------------
void spectrum_array::init(const std::vector<double> w,
		const int n_mu, const int n_phi)
{
	// assign wave grid
	double w_start = w[0];
	double w_stop  = w[1];
	double w_del   = w[2];
	nu_grid.init(w_start,w_stop,w_del);
	int n_wave   = nu_grid.size();

	// asign mu grid
	mu_grid.init(-1,1,n_mu);

	// asign phi grid
	phi_grid.init(-pc::pi,pc::pi,n_phi);

	// index parameters
	unsigned n_elements  = n_wave*n_mu*n_phi;

	// allocate memory
	flux.resize(n_elements);

	// clear
	wipe();
}


//--------------------------------------------------------------
// Initialization and Allocation
//--------------------------------------------------------------
void spectrum_array::init(const locate_array wg,
		const locate_array mg, const locate_array pg)
{
	// initialize locate arrays by swapping with the inputs
	nu_grid.copy(wg);
	mu_grid.copy(mg);
	phi_grid.copy(pg);

	int n_wave   = nu_grid.size();
	int n_mu     = mu_grid.size();
	int n_phi    = phi_grid.size();

	// index parameters
	unsigned n_elements  = n_wave*n_mu*n_phi;

	// allocate memory
	flux.resize(n_elements);

	// clear
	wipe();
}


//--------------------------------------------------------------
// Functional procedure: Wipe
//--------------------------------------------------------------
void spectrum_array::wipe()
{
	for(unsigned i=0;i<flux.size();i++) flux[i] = 0;
}



//--------------------------------------------------------------
// handles the indexing: should be called in this order
//    group, mu, phi
//--------------------------------------------------------------
unsigned spectrum_array::index(const unsigned nu_bin, const unsigned mu_bin, const unsigned phi_bin) const
{
	unsigned n_phi = phi_grid.size();
	unsigned n_nu  =  nu_grid.size();
	unsigned n_mu  =  mu_grid.size();
	PRINT_ASSERT(nu_bin,<,n_nu);
	PRINT_ASSERT(mu_bin,<,n_mu);
	PRINT_ASSERT(phi_bin,<,n_phi);
	int a3 = n_phi;
	int a2 = n_mu*a3;
	const unsigned ind = nu_bin*a2 + mu_bin*a3 + phi_bin;
	PRINT_ASSERT(ind,<,n_phi*n_nu*n_mu);
	return ind;
}

//----------------
// get bin indices
//----------------
unsigned spectrum_array::nu_bin(const unsigned index) const{
	PRINT_ASSERT(index,<,flux.size());
	unsigned result = index / (phi_grid.size()*mu_grid.size());
	PRINT_ASSERT(result,<,nu_grid.size());
	return result;
}
unsigned spectrum_array::mu_bin(const unsigned index) const{
	PRINT_ASSERT(index,<,flux.size());
	double result =  (index%(phi_grid.size()*mu_grid.size())) / phi_grid.size();
	PRINT_ASSERT(result,<,mu_grid.size());
	return result;
}
unsigned spectrum_array::phi_bin(const unsigned index) const{
	PRINT_ASSERT(index,<,flux.size());
	double result = index % phi_grid.size();
	PRINT_ASSERT(result,<,phi_grid.size());
	return result;
}

//--------------------
// get centers of bins
//--------------------
double spectrum_array::nu_center(const unsigned index) const{
	PRINT_ASSERT(index,<,flux.size());
	return nu_grid.center(nu_bin(index));
}
double spectrum_array::mu_center(const unsigned index) const{
	PRINT_ASSERT(index,<,flux.size());
	return mu_grid.center(mu_bin(index));
}
double spectrum_array::phi_center(const unsigned index) const{
	PRINT_ASSERT(index,<,flux.size());
	return phi_grid.center(phi_bin(index));
}
double spectrum_array::nu_bin_center(const unsigned index) const{
	PRINT_ASSERT(index,<,nu_grid.size());
	return nu_grid.center(index);
}
double spectrum_array::mu_bin_center(const unsigned index) const{
	PRINT_ASSERT(index,<,mu_grid.size());
	return mu_grid.center(index);
}
double spectrum_array::phi_bin_center(const unsigned index) const{
	PRINT_ASSERT(index,<,phi_grid.size());
	return phi_grid.center(index);
}

//-------
// getter
//-------
double spectrum_array::get(const int index) const{
	return flux[index];
}

//----------------------
// get size of spectrum
//----------------------
unsigned spectrum_array::size() const{
	return flux.size();
}

//--------------------------------------------------------------
// count a particle
////--------------------------------------------------------------
void spectrum_array::count(const double D[3], const int Dsize, const double nu, const double E)
{
	PRINT_ASSERT(Dsize,==,3);
	PRINT_ASSERT(E,>=,0);
	const double tiny = 1e-8;
	double mu  = D[2];           // component along z axis
	mu = max(-1.0+tiny,mu);
	mu = min( 1.0-tiny,mu);
	double phi = atan2(D[1],D[0]);  // projection into x-y plane
	if(phi< -pc::pi) phi += 2.0*pc::pi;
	if(phi>  pc::pi) phi -= 2.0*pc::pi;

	// if off the LEFT of mu/phi grids, just return without counting
	if ((mu<mu_grid.min) || (phi<phi_grid.min)){
		cout << "Lost a particle off the left of the spectrum array!" << endl;
		return;
	}

	// locate bin number in all dimensions.
	unsigned nu_bin  =  nu_grid.locate(nu);
	unsigned mu_bin  =  mu_grid.locate(mu);
	unsigned phi_bin = phi_grid.locate(phi);

	// if off the RIGHT of mu/phi grids, just return without counting
	if((mu_bin ==   mu_grid.size()) || (phi_bin ==  phi_grid.size())){
		cout << "Lost a particle off the right of the spectrum array!" << endl;
		return;
	}

	// if off RIGHT of wavelength grid, store in last bin (LEFT is accounted for by locate)
	if (nu_bin == nu_grid.size()) nu_bin--;

	// add to counters
	int ind = index(nu_bin,mu_bin,phi_bin);

    #pragma omp atomic
	flux[ind]  += E;
	PRINT_ASSERT(flux[ind],>=,0);
	PRINT_ASSERT(E,<,INFINITY);
	PRINT_ASSERT(flux[ind],<,INFINITY);
}

double spectrum_array::average_nu() const{
	double integral1 = 0;
	double integral2 = 0;
	for(unsigned nu_bin=0; nu_bin<nu_grid.size(); nu_bin++){
		for(unsigned mu_bin=0; mu_bin<mu_grid.size(); mu_bin++){
			for(unsigned phi_bin=0; phi_bin<phi_grid.size(); phi_bin++){
				int ind = index(nu_bin,mu_bin,phi_bin);
				integral1 += flux[ind];
				integral2 += flux[ind] * nu_grid.center(nu_bin);
			}
		}
	}
	if(integral2/integral1 >= 0) return integral2 / integral1;
	else return 0;
}

double spectrum_array::integrate() const{
	double integral = 0;
	for(unsigned nu_bin=0; nu_bin<nu_grid.size(); nu_bin++){
		for(unsigned mu_bin=0; mu_bin<mu_grid.size(); mu_bin++){
			for(unsigned phi_bin=0; phi_bin<phi_grid.size(); phi_bin++){
				int ind = index(nu_bin,mu_bin,phi_bin);
				integral += flux[ind];
			}
		}
	}
	return integral;
}

// integrate over direction
void spectrum_array::integrate_over_direction(vector<double>& integral) const{
	integral = vector<double>(nu_grid.size(),0);
	for(unsigned nu_bin=0; nu_bin<nu_grid.size(); nu_bin++){
		for(unsigned mu_bin=0; mu_bin<mu_grid.size(); mu_bin++){
			for(unsigned phi_bin=0; phi_bin<phi_grid.size(); phi_bin++){
				int ind = index(nu_bin,mu_bin,phi_bin);
				integral[nu_bin] += flux[ind];
			}
		}
	}
}

//--------------------------------------------------------------
// print out
//--------------------------------------------------------------
void spectrum_array::print(const int iw, const int species) const
{
	ofstream outf;
	stringstream speciesstream;
	speciesstream << species;
	string prefix = "spectrum_species"+speciesstream.str();
	string filename = transport::filename(prefix.c_str(),iw,".dat");
	outf.open(filename.c_str());

	unsigned n_nu  =  nu_grid.size();
	unsigned n_mu  =  mu_grid.size();
	unsigned n_phi = phi_grid.size();
	double solid_angle = 4.0*pc::pi / (double)(n_mu*n_phi);

	outf << "# " << "n_nu:" << n_nu << " " << "n_mu:" << n_mu << " " << "n_phi:" << n_phi << endl;
	outf << "# ";
	outf << (n_nu >1 ? "frequency(Hz) delta(Hz)" : "");
	outf << (n_mu >1 ? "mu "            : "");
	outf << (n_phi>1 ? "phi "           : "");
	outf << "integrated_flux(erg/s/Hz/sr)" << endl;

	for (unsigned k=0;k<n_mu;k++)
		for (unsigned m=0;m<n_phi;m++)
			for (unsigned j=0;j<n_nu;j++)
			{
				int id = index(j,k,m);
				if (n_nu > 1) outf <<  nu_grid.center(j) << " " << nu_grid.delta(j) << " ";
				if (n_mu > 1) outf <<  mu_grid.center(k) << " ";
				if (n_phi> 1) outf << phi_grid.center(m) << " ";

				// the delta is infinity if the bin is a catch-all.
				double wdel = (nu_grid.delta(j)<numeric_limits<double>::infinity() ? nu_grid.delta(j) : 1);
				outf << flux[id]/wdel/solid_angle << endl;
			}
	outf.close();
}


void  spectrum_array::rescale(double r)
{
	for(unsigned i=0;i<flux.size();i++) flux[i] *= r;
}


//--------------------------------------------------------------
// MPI average the spectrum contents
//--------------------------------------------------------------
// only process 0 gets the reduced spectrum to print
void spectrum_array::MPI_average()
{
	const unsigned n_elements = nu_grid.size()*mu_grid.size()*phi_grid.size();

	// average the flux (receive goes out of scope after section)
	vector<double> receive;
	receive.resize(n_elements);
	MPI_Allreduce(&flux.front(), &receive.front(), n_elements, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	for(unsigned i=0; i<flux.size(); i++) flux[i] = receive[i];//flux.swap(receive);

	// only have the receiving ID do the division
	int myID, mpi_procs;
	MPI_Comm_size( MPI_COMM_WORLD, &mpi_procs );
	MPI_Comm_rank( MPI_COMM_WORLD, &myID      );
	rescale(1./(double)mpi_procs);
}

//--------------------------------------------------------------
// Write data to specified location in an HDF5 file
//--------------------------------------------------------------
void spectrum_array::write_hdf5_data(H5::DataSet dataset, H5::DataSpace dataspace) const
{
	// some sanity checks
	unsigned ndims = dataspace.getSimpleExtentNdims();
	vector<hsize_t> dspace_dims(ndims,0.0);
	dataspace.getSimpleExtentDims(&dspace_dims[0]);
	PRINT_ASSERT(ndims,>=,4);
	PRINT_ASSERT(dspace_dims[ndims-4],<=,6);
	PRINT_ASSERT(dspace_dims[ndims-3],==,nu_grid.size());
	PRINT_ASSERT(dspace_dims[ndims-2],==,mu_grid.size());
	PRINT_ASSERT(dspace_dims[ndims-1],==,phi_grid.size());

	// define the memory dataspace
	hsize_t mdim[3];
	mdim[0] =  nu_grid.size();
	mdim[1] =  mu_grid.size();
	mdim[2] = phi_grid.size();
	H5::DataSpace memspace(3,mdim,NULL);

	// write the data (converting to single precision)
	// assumes phi increases fastest, then mu, then nu
	vector<float> tmp(flux.size(),-1.0);
	for(unsigned i=0; i<flux.size(); i++) tmp[i] = flux[i];
	dataset.write(&tmp[0], H5::PredType::IEEE_F32LE, memspace, dataspace);
}

//--------------------------------------------------------------
// Write distribution function coordinates to an HDF5 file
//--------------------------------------------------------------
void spectrum_array::write_hdf5_coordinates(H5::H5File file) const
{
	// useful quantities
	hsize_t dims[1];
	H5::DataSpace dataspace;
	H5::DataSet dataset;
	vector<float> tmp;

	// write nu_grid
	tmp = vector<float>(nu_grid.size()+1,0.0);
	dims[0] = nu_grid.size()+1;
	dataspace = H5::DataSpace(1,dims);
	dataset = file.createDataSet("distribution_frequency_grid(Hz,lab)",H5::PredType::IEEE_F32LE,dataspace);
	tmp[0] = nu_grid.min;
	for(unsigned i=1; i<nu_grid.size()+1; i++) tmp[i] = nu_grid[i-1];
	dataset.write(&tmp[0],H5::PredType::IEEE_F32LE);
	dataset.close();

	// write mu_grid
	tmp = vector<float>(mu_grid.size()+1,0.0);
	dims[0] = mu_grid.size()+1;
	dataspace = H5::DataSpace(1,dims);
	dataset = file.createDataSet("distribution_costheta_grid(lab)",H5::PredType::IEEE_F32LE,dataspace);
	tmp[0] = mu_grid.min;
	for(unsigned i=1; i<mu_grid.size()+1; i++) tmp[i] = mu_grid[i-1];
	dataset.write(&tmp[0],H5::PredType::IEEE_F32LE);
	dataset.close();

	// write phi_grid
	tmp = vector<float>(phi_grid.size()+1,0.0);
	dims[0] = phi_grid.size()+1;
	dataspace = H5::DataSpace(1,dims);
	dataset = file.createDataSet("distribution_phi_grid(radians,lab)",H5::PredType::IEEE_F32LE,dataspace);
	tmp[0] = phi_grid.min;
	for(unsigned i=1; i<phi_grid.size()+1; i++) tmp[i] = phi_grid[i-1];
	dataset.write(&tmp[0],H5::PredType::IEEE_F32LE);
	dataset.close();
}
