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

#ifndef _MOMENT_SPECTRUM_ARRAY_H
#define _MOMENT_SPECTRUM_ARRAY_H 1

#include "SpectrumArray.h"
#include "H5Cpp.h"
#include <fstream>

using namespace std;

template<unsigned ndims_spatial>
class MomentSpectrumArray : public SpectrumArray {

private:

	// helper variables
	static const unsigned index_range = 3; // index can be 0,1,2
	static const unsigned nranks = 4;
	static const unsigned nuGridIndex = ndims_spatial;
	const unsigned n_independent_elements[4] = {1,3,6,10}; // for a symmetric tensor of rank 0,1,2,3
	static const unsigned n_total_elements = 1+3+6+10;

	// bin arrays
	// values represent bin upper walls (the single locate_array.min value is the leftmost wall)
	// underflow is combined into leftmost bin (right of the locate_array.min)
	// overflow is combined into the rightmost bin (left of locate_array[size-1])
	MultiDArray<n_total_elements, ndims_spatial+1> data;

public:

	//---------------------------------------------------
	// increment tensor indices for any symmetric tensor
	//---------------------------------------------------
	static void increment_tensor_indices(unsigned tensor_indices[], const unsigned rank) {
		PRINT_ASSERT(rank,>=,0);
		PRINT_ASSERT(rank,<,nranks);
		if (rank == 0)
			return;
		tensor_indices[0]++;
		for (int which_index = 0; which_index < rank - 1; which_index++) {
			if (tensor_indices[which_index] > index_range - 1) {
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

	void init(const vector<Axis>& spatial_axes, const Axis& nu_grid) {
		PRINT_ASSERT(spatial_axes.size(),==,ndims_spatial);
		vector<Axis> axes(ndims_spatial+1);

		// spatial axes
		for(int i=0; i<ndims_spatial; i++) axes[i] = spatial_axes[i];

		// frequency axis
		axes[nuGridIndex] = nu_grid;

		// initialize the moments
		data.set_axes(axes);
		wipe();
	}

	//--------------------------------------------------------------
	// Functional procedure: Wipe
	//--------------------------------------------------------------
	void wipe(){
		data.wipe();
	}

	//--------------------------------------------------------------
	// count a particle
	////--------------------------------------------------------------
	void count(const EinsteinHelper* eh, const double E) {
		PRINT_ASSERT(E, >=, 0);
		PRINT_ASSERT(E, <, INFINITY);
		double D[3] = {eh->kup_tet[0], eh->kup_tet[1], eh->kup_tet[2]};
		Metric::normalize_Minkowski<3>(D);
		
		unsigned data_indices[data.Ndims()];
		for(int i=0; i<ndims_spatial; i++) data_indices[i] = eh->dir_ind[i];
		data_indices[nuGridIndex] = eh->dir_ind[NDIMS];

		// increment moments
		Tuple<double, n_total_elements> tmp;
		unsigned tuple_index=0;
		for(unsigned rank = 0; rank<nranks; rank++) {
			unsigned tensor_indices[rank];
			for(unsigned r = 0; r<rank; r++) tensor_indices[r] = 0;
			for(unsigned i=0; i<n_independent_elements[rank]; i++) {
				tmp[tuple_index] = E;
				for(int r = 0; r<rank; r++) tmp[tuple_index] *= D[tensor_indices[r]];
				increment_tensor_indices(tensor_indices, rank);
				tuple_index++;
			}
		}
		data.add(data_indices, tmp);
	}


	void rescale(const double r) {
		for(unsigned i=0; i<data.size(); i++) data[i] *= r;
	}
	void rescale_spatial_point(const unsigned dir_ind[ndims_spatial], const double r){
		unsigned all_indices[ndims_spatial+1];
		for(unsigned i=0; i<ndims_spatial; i++) all_indices[i] = dir_ind[i];
		all_indices[ndims_spatial  ] = 0;
		unsigned base_ind = data.direct_index(all_indices);
		unsigned nbins = data.axes[ndims_spatial].size();
		for(unsigned i=0; i<nbins; i++){
			data.y0[base_ind+i] *= r;
		}
	}

	//--------------------------------------------------------------
	// MPI average the spectrum contents
	//--------------------------------------------------------------
	// only process 0 gets the reduced spectrum to print
	void MPI_average() {
		data.MPI_combine();
	}

	//--------------------------------------------------------------
	// Write data to specified location in an HDF5 file
	//--------------------------------------------------------------
	void write_hdf5_data(H5::H5File file,const string name) const {
		data.write_HDF5(file, name);
	}

	//--------------------------------------------------------------
	// Write distribution function coordinates to an HDF5 file
	//--------------------------------------------------------------
	void write_hdf5_coordinates(H5::H5File file,
			const string name) const {
		// useful quantities
		hsize_t dims[1];
		H5::DataSpace dataspace;
		H5::DataSet dataset;
		H5::Group group;
		vector<float> tmp;

		// set up the dfunc group
		stringstream datasetname, indicesname;

		// SET UP DATASPACE FOR EACH MOMENT
		for (unsigned rank=0; rank<nranks; rank++) {
			// prep the filenames
			indicesname.str("");
			indicesname << name;
			indicesname << "_moment_" << rank << "_indices";

			// set up the database with the indices
			hsize_t indices_dimensions[] = { n_independent_elements[rank], rank };
			dataspace = H5::DataSpace(2, indices_dimensions);
			dataset = file.createDataSet(indicesname.str(), H5::PredType::STD_I32LE,
					dataspace);
			unsigned tensor_indices[rank];
			for (int r = 0; r < rank; r++)
				tensor_indices[r] = 0;
			int indices2D[n_independent_elements[rank]][rank];
			for (int i = 0; i < n_independent_elements[rank]; i++) {
				for (int r = 0; r < rank; r++)
					indices2D[i][r] = tensor_indices[r];
				increment_tensor_indices(tensor_indices, rank);
			}
			dataset.write(indices2D, H5::PredType::STD_I32LE);
			dataset.close();

		}
	}
	void add_isotropic(const unsigned dir_ind[NDIMS+1], const double E){
		unsigned indices[data.Ndims()];
		for(int i=0; i<ndims_spatial; i++) indices[i] = dir_ind[i];
		indices[nuGridIndex] = dir_ind[ndims_spatial];

		// increment moments
		Tuple<double, n_total_elements> tmp;
		tmp = 0;
		tmp[0] = E;
		tmp[4] = E/3.; // xx
		tmp[7] = E/3.; // yy
		tmp[9] = E/3.; // zz
		data.add(indices, tmp);
	}

	double total() const{
		double result=0;
		for(unsigned i=0; i<data.size(); i++)
			result += data[i][0];
		return result;
	}

};

#endif
