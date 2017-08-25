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
#include <sstream>
#include <fstream>
#include "global_options.h"
#include "RadialMomentSpectrumArray.h"
#include "Transport.h"

using namespace std;
namespace pc = physical_constants;


//--------------------------------------------------------------
// Initialization and Allocation
//--------------------------------------------------------------
void RadialMomentSpectrumArray::init(const LocateArray wg, const int max_rank) {
	// initialize locate arrays by swapping with the inputs
	nu_grid.copy(wg);

	// initialize the moments
	nranks = max_rank + 1;
	moments.resize(nu_grid.size() * nranks);
	
	// clear the data
	wipe();
}

//--------------------------------------------------------------
// Functional procedure: Wipe
//--------------------------------------------------------------
void RadialMomentSpectrumArray::wipe() {
	for (int i = 0; i < moments.size(); i++)
		moments[i] = 0;
}

//----------------------
// get bin central value
//----------------------
double RadialMomentSpectrumArray::nu_bin_center(const unsigned index) const {
	PRINT_ASSERT(index, <, nu_grid.size());
	return nu_grid.center(index);
}

//----------------------
// get size of spectrum
//----------------------
unsigned RadialMomentSpectrumArray::index(const int nu_bin, const int rank) const {
	return nranks*nu_bin + rank;
}


//--------------------------------------------------------------
// count a particle
////--------------------------------------------------------------
void RadialMomentSpectrumArray::count(const double D[3], const int Dsize,
		const double nu, const double E) {
	PRINT_ASSERT(Dsize, ==, 3);
	PRINT_ASSERT(E, >=, 0);
	PRINT_ASSERT(E, <, INFINITY);

	// locate bin number
	unsigned nu_bin = nu_grid.locate(nu);
	if (nu_bin == nu_grid.size()) nu_bin--;

	// increment moments
	double tmp;
	for (int rank = 0; rank<nranks; rank++) {
		int index2D = index(nu_bin,rank);
		tmp = E;
		for (int r = 0; r<rank; r++) tmp *= D[2];
		#pragma omp atomic
		moments[index2D] += tmp;
	}
}

double RadialMomentSpectrumArray::average_nu() const {
	double integral_E = 0;
	double integral_N = 0;
	for (unsigned nu_bin = 0; nu_bin<nu_grid.size(); nu_bin++) {
		int index2D = index(nu_bin,0);
		integral_E += moments[index2D];
		integral_N += moments[index2D] / nu_grid.center(nu_bin);
	}
	double result = integral_E / integral_N;
	if (result >= 0)
		return result;
	else
		return 0;
}

double RadialMomentSpectrumArray::integrate() const {
	double integral = 0;
	for (unsigned nu_bin = 0; nu_bin<nu_grid.size(); nu_bin++){
		int index2D = index(nu_bin,0);
		integral += moments[index2D];
	}
	return integral;
}

// integrate over direction
void RadialMomentSpectrumArray::integrate_over_direction(
		vector<double>& integral) const {
	integral = vector<double>(nu_grid.size(), 0);
	for (unsigned nu_bin = 0; nu_bin<nu_grid.size(); nu_bin++) {
		int index2D = index(nu_bin,0);
		integral[nu_bin] = moments[index2D];
	}
}

//--------------------------------------------------------------
// print out
//--------------------------------------------------------------
void RadialMomentSpectrumArray::print(const int iw, const int species) const {
	ofstream outf;
	stringstream speciesstream;
	speciesstream << species;
	string prefix = "spectrum_species" + speciesstream.str();
	string filename = Transport::filename(prefix.c_str(), iw, ".dat");
	outf.open(filename.c_str());

	// write the header
	outf << "# " << "n_nu:" << nu_grid.size() << " n_ranks:" << nranks << endl;
	outf << "# ";
	outf << "frequency(Hz) delta(Hz) M0(erg/s/Hz)";
	for (int rank = 1; rank < nranks; rank++) {
		int tensor_indices[rank];
		for (int r = 0; r < rank; r++)
			tensor_indices[r] = 0;
		while (tensor_indices[rank - 1] < 3) {
			outf << " M" << rank << ":";
			for (int r = 0; r < rank; r++)
				outf << tensor_indices[r];
		}
	}
	outf << endl;

	// write the data
	for (unsigned group = 0; group < nu_grid.size(); group++) {
		outf << nu_grid.center(group) << " " << nu_grid.delta(group) << " ";
		for (unsigned rank = 0; rank < nranks; rank++){
			int index2D = index(group, rank);
			// the delta is infinity if the bin is a catch-all.
			double wdel = (nu_grid.delta(group)<numeric_limits<double>::infinity() ? nu_grid.delta(group) : 1);
			outf << moments[index2D] / wdel << endl;
		}
	}
	outf.close();
}

void RadialMomentSpectrumArray::rescale(double r) {
	for(int index2D=0; index2D<moments.size(); index2D++)
		moments[index2D] *= r;
}


//--------------------------------------------------------------
// MPI average the spectrum contents
//--------------------------------------------------------------
// only process 0 gets the reduced spectrum to print
void RadialMomentSpectrumArray::MPI_average(const int proc)
{
	int MPI_myID;
	MPI_Comm_rank( MPI_COMM_WORLD, &MPI_myID      );
	MPI_Request request;

	const unsigned n_elements = moments.size();
	const int tag=0;
	
	// average the flux (receive goes out of scope after section)
	vector<double> receive;
	receive.resize(n_elements);
	MPI_Reduce(&moments.front(), &receive.front(), n_elements, MPI_DOUBLE, MPI_SUM, proc, MPI_COMM_WORLD);
	if(proc>0){
	  if(MPI_myID==0) MPI_Irecv(&receive.front(), n_elements, MPI_DOUBLE, proc, tag, MPI_COMM_WORLD, &request);
	  if(MPI_myID==proc) MPI_Isend(&receive.front(), n_elements, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &request);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	for(unsigned i=0; i<moments.size(); i++) moments[i] = receive[i];//flux.swap(receive);
}


//--------------------------------------------------------------
// Write data to specified location in an HDF5 file
//--------------------------------------------------------------
void RadialMomentSpectrumArray::write_hdf5_data(H5::H5File file, const int s,
		const int dir_ind[], const hsize_t n_spatial_dims) const {

	// get the dataset
	H5::DataSet dataset = file.openDataSet("distribution(erg|ccm,lab)");

	// get the dataspace of the dataset
	H5::DataSpace dataspace = dataset.getSpace();
	hsize_t n_total_dims = dataspace.getSimpleExtentNdims();
	hsize_t dfunc_dims[n_spatial_dims+3];
	hsize_t start[n_spatial_dims+3];
	dataspace.getSelectBounds(start,dfunc_dims);
	for(int i=0; i<n_total_dims; i++) dfunc_dims[i] += 1;
	PRINT_ASSERT(dataspace.getSimpleExtentNdims(),==,n_spatial_dims+3);
	PRINT_ASSERT(dfunc_dims[n_spatial_dims+1],==,nu_grid.size());
	PRINT_ASSERT(dfunc_dims[n_spatial_dims+2],==,nranks);

	// set up subspace offset
	hsize_t offset[n_total_dims];
	for(unsigned i=0; i<n_spatial_dims; i++) offset[i] = dir_ind[i];
	offset[n_spatial_dims  ] = s;
	offset[n_spatial_dims+1] = 0;
	offset[n_spatial_dims+2] = 0;

	// set up subspace stride
	hsize_t stride[n_total_dims];
	for(unsigned i=0; i<n_total_dims; i++) stride[i] = 1;

	// set up subspace block
	hsize_t block[n_total_dims];
	for(unsigned i=0; i<n_total_dims; i++) block[i] = 1;

	// set up local spectrum dimensions
	hsize_t spectrum_dims[n_total_dims];
	for(unsigned i=0; i<n_spatial_dims; i++) spectrum_dims[i] = 1;
	spectrum_dims[n_spatial_dims  ] = 1;
	spectrum_dims[n_spatial_dims+1] = nu_grid.size();
	spectrum_dims[n_spatial_dims+2] = nranks;

	// set dataspace
	dataspace.selectHyperslab(H5S_SELECT_SET,&spectrum_dims[0],&offset[0],&stride[0],&block[0]);

	// define the memory dataspace
	hsize_t mdim[2];
	mdim[0] =  nu_grid.size();
	mdim[1] =  nranks;
	H5::DataSpace memspace(2,mdim,NULL);

	// write the data (converting to single precision)
	// assumes phi increases fastest, then mu, then nu
	vector<float> tmp(moments.size(),-1.0);
	for(unsigned i=0; i<moments.size(); i++) tmp[i] = moments[i];
	dataset.write(&tmp[0], H5::PredType::IEEE_F32LE, memspace, dataspace);
	dataset.close();
}

//--------------------------------------------------------------
// Write distribution function coordinates to an HDF5 file
//--------------------------------------------------------------
void RadialMomentSpectrumArray::write_hdf5_coordinates(H5::H5File file,
		const Grid* grid) const {

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

	// SET UP +3D DATASPACE
	hsize_t zdims[grid->dimensionality()];
	grid->dims(zdims,grid->dimensionality());
	vector<hsize_t> dims_plus3(grid->dimensionality()+3,0);
	for(unsigned i=0; i<grid->dimensionality(); i++) dims_plus3[i] = zdims[i]; // number of spatial bins
	dims_plus3[grid->dimensionality()  ] = grid->z[0].distribution.size(); // number of species
	dims_plus3[grid->dimensionality()+1] = nu_grid.size(); // number of energy bins
	dims_plus3[grid->dimensionality()+2] = nranks; // number of moments
	dataspace = H5::DataSpace(grid->dimensionality()+3,&dims_plus3[0]);
	dataset = file.createDataSet("distribution(erg|ccm,lab)",H5::PredType::IEEE_F32LE,dataspace);
}

void RadialMomentSpectrumArray::write_header(ofstream& outf) const {
	for (int group = 0; group < nu_grid.size(); group++)
		for (int rank = 0; rank < nranks; rank++)
			outf << "ig" << group << "M" << rank << "(erg/ccm)  ";
}

void RadialMomentSpectrumArray::write_line(ofstream& outf) const {
	for (int index2D = 0; index2D < moments.size(); index2D++)
		outf << moments[index2D] << "\t";
}
