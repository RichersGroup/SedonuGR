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
#include "MomentSpectrumArray.h"
#include "Transport.h"

using namespace std;
namespace pc = physical_constants;

//---------------------------------------------------
// increment tensor indices for any symmetric tensor
//---------------------------------------------------
void increment_tensor_indices(int tensor_indices[], const int rank, const int length) {
	PRINT_ASSERT(rank, >=, 0);
	PRINT_ASSERT(length, >, 0);
	if (rank == 0)
		return;
	tensor_indices[0]++;
	for (int which_index = 0; which_index < rank - 1; which_index++) {
		if (tensor_indices[which_index] > length - 1) {
			tensor_indices[which_index + 1]++;
			for (int which_index_2 = 0; which_index_2 < which_index + 1;
					which_index_2++)
				tensor_indices[which_index_2] = tensor_indices[which_index + 1];
		}
	}
}

//--------------------------------------------------------------
// Initialization and Allocation
//--------------------------------------------------------------
void MomentSpectrumArray::init(const vector<Axis>& spatial_axes, const Axis& nu_grid, const int max_rank) {
	vector<Axis> axes(spatial_axes.size()+2);
	nuGridIndex = spatial_axes.size()+1;
	momGridIndex = spatial_axes.size()+2;

	// spatial axes
	for(int i=0; i<spatial_axes.size(); i++) axes[i] = spatial_axes[i];

	// frequency axis
	axes[nuGridIndex] = nu_grid;


	// determine the number of independent elements
	const int length = 3;
	int n_independent_elements[max_rank + 1];
	n_independent_elements[0] = 1;
	for (int rank = 1; rank < max_rank + 1; rank++) {
		n_independent_elements[rank] = 0;
		int tensor_indices[rank];
		for (int i = 0; i < rank; i++)
			tensor_indices[i] = 0;
		while (tensor_indices[rank - 1] < length) {
			n_independent_elements[rank]++;
			increment_tensor_indices(tensor_indices, rank, length);
		}
	}

	// initialize the moments
	nranks = max_rank+1;
	data.resize(nranks);
	for (int rank=0; rank<nranks; rank++) {
		double minval = -0.5; // this whole thing is a hack. do not try to use to interpolate...
		unsigned nmoments = n_independent_elements[rank];
		vector<double> top(nmoments);
		vector<double> mid(nmoments);
		for(int i=0; i<nmoments; i++){
			mid[i] = i;
			top[i] = i + 0.5;
		}
		axes[momGridIndex] = Axis(minval, top, mid);

		switch(spatial_axes.size()){
		case 0: data[rank] = new MultiDArray<0+2>(axes); break;
		case 1: data[rank] = new MultiDArray<1+2>(axes); break;
		case 2: data[rank] = new MultiDArray<2+2>(axes); break;
		case 3: data[rank] = new MultiDArray<3+2>(axes); break;
		}

		data[rank]->wipe();
	}
}

//--------------------------------------------------------------
// Functional procedure: Wipe
//--------------------------------------------------------------
void MomentSpectrumArray::wipe() {
	for(int rank=0; rank<nranks; rank++) data[rank]->wipe();
}

//--------------------------------------------------------------
// count a particle
////--------------------------------------------------------------
void MomentSpectrumArray::count(const double D[3], const vector<unsigned>& dir_ind, const double nu, const double E) {
	PRINT_ASSERT(E, >=, 0);
	PRINT_ASSERT(E, <, INFINITY);

	unsigned data_indices[data[0]->Ndims()];
	for(int i=0; i<dir_ind.size(); i++) data_indices[i] = dir_ind[i];

	int nu_bin = data[0]->axes[nuGridIndex].bin(nu);
	nu_bin = max(nu_bin, 0);
	nu_bin = min(nu_bin, data[0]->axes[nuGridIndex].size()-1);
	data_indices[nuGridIndex] = nu_bin;

	// increment moments
	double tmp;
	for (int rank = 0; rank<nranks; rank++) {
		const int nelements = data[rank]->axes[momGridIndex].size();
		int tensor_indices[rank];
		for (int r = 0; r<rank; r++) tensor_indices[r] = 0;
		
		for (int i = 0; i<nelements; i++) {
			data_indices[momGridIndex] = i;
			tmp = E;
			for (int r = 0; r<rank; r++) tmp *= D[tensor_indices[r]];
			increment_tensor_indices(tensor_indices, rank, 3);

			data[rank]->add(data_indices, tmp);
		}
	}
}



void MomentSpectrumArray::rescale(double r) {
	for(unsigned rank=0; rank<data.size(); rank++)
		for(unsigned i=0; i<data[i]->size(); i++)
			data[rank]->y0[i] *= r;
}

//--------------------------------------------------------------
// MPI average the spectrum contents
//--------------------------------------------------------------
// only process 0 gets the reduced spectrum to print
void MomentSpectrumArray::MPI_average() {
	for(unsigned rank=0; rank<data.size(); rank++)
		data[rank]->MPI_combine();
}

//--------------------------------------------------------------
// Write data to specified location in an HDF5 file
//--------------------------------------------------------------
void MomentSpectrumArray::write_hdf5_data(H5::H5File file,const string name) const {
	for (int rank = 0; rank<data.size(); rank++)
		data[rank]->write_HDF5(file, name+"[r"+to_string(rank)+"]");
}

//--------------------------------------------------------------
// Write distribution function coordinates to an HDF5 file
//--------------------------------------------------------------
void MomentSpectrumArray::write_hdf5_coordinates(H5::H5File file,
		const string name) const {
	// useful quantities
	hsize_t dims[1];
	H5::DataSpace dataspace;
	H5::DataSet dataset;
	H5::Group group;
	vector<float> tmp;

	// set up the dfunc group
	string dfunc_dir = "distribution(erg|ccm,lab)";
	group = file.createGroup(dfunc_dir);
	stringstream datasetname, indicesname;

	// SET UP DATASPACE FOR EACH MOMENT
	for (unsigned rank = 0; rank < data.size(); rank++) {
		unsigned n_elements = data[rank]->axes[momGridIndex].size();
		// prep the filenames
		indicesname.str("");
		indicesname << "moment_" << rank << "_indices";

		// set up the database with the indices
		hsize_t indices_dimensions[] = { n_elements, rank };
		dataspace = H5::DataSpace(2, indices_dimensions);
		dataset = file.createDataSet(indicesname.str(), H5::PredType::STD_I32LE,
				dataspace);
		int tensor_indices[rank];
		for (int r = 0; r < rank; r++)
			tensor_indices[r] = 0;
		int indices2D[n_elements][rank];
		for (int i = 0; i < n_elements; i++) {
			for (int r = 0; r < rank; r++)
				indices2D[i][r] = tensor_indices[r];
			increment_tensor_indices(tensor_indices, rank, 3);
		}
		dataset.write(indices2D, H5::PredType::STD_I32LE);
		dataset.close();

	}
}

