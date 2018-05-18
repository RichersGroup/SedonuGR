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

#ifndef _RADIAL_MOMENT_SPECTRUM_ARRAY_H
#define _RADIAL_MOMENT_SPECTRUM_ARRAY_H 1

#include "SpectrumArray.h"
#include "H5Cpp.h"
#include <fstream>

using namespace std;

template<unsigned ndims_spatial>
class RadialMomentSpectrumArray : public SpectrumArray {

private:

	// bin arrays
	// values represent bin upper walls (the single locate_array.min value is the leftmost wall)
	// underflow is combined into leftmost bin (right of the locate_array.min)
	// overflow is combined into the rightmost bin (left of locate_array[size-1])
	MultiDArray<double,4,ndims_spatial+1> data; // 0, r, rr, rrr

	static const unsigned nranks = 4;
	static const unsigned nuGridIndex = ndims_spatial;

public:

	//--------------------------------------------------------------
	// Initialization and Allocation
	//--------------------------------------------------------------
	void init(const vector<Axis>& spatial_axes, const Axis& nu_grid) {
		vector<Axis> axes;
		for(unsigned i=0; i<spatial_axes.size(); i++) axes.push_back(spatial_axes[i]);

		axes.push_back(nu_grid);

		// set up the data structure
		data.set_axes(axes);
		data.wipe();
	}

	//--------------------------------------------------------------
	// Functional procedure: Wipe
	//--------------------------------------------------------------
	void wipe() {
		data.wipe();
	}

	//--------------------------------------------------------------
	// count a particle
	////--------------------------------------------------------------
	void count(const Tuple<double,4>& kup_tet, const unsigned dir_ind[NDIMS+1], const double E) {
		PRINT_ASSERT(E, >=, 0);
		PRINT_ASSERT(E, <, INFINITY);
		Tuple<double,3> D = {kup_tet[0], kup_tet[1], kup_tet[2]};
		Metric::normalize_Minkowski<3>(D);

		unsigned indices[data.Ndims()];
		for(unsigned i=0; i<ndims_spatial; i++) indices[i] = dir_ind[i];
		indices[nuGridIndex] = dir_ind[NDIMS];

		// increment moments
		Tuple<double,4> tmp;
		tmp[0] = E;
		for (unsigned rank = 1; rank<nranks; rank++) {
			tmp[rank] = tmp[rank-1]*D[2];
		}
		data.add(indices, tmp);
	}

	void rescale(double r) {
		for(unsigned i=0;i<data.size();i++) data.y0[i] *= r;
	}
	void rescale_spatial_point(const unsigned dir_ind[ndims_spatial], const double r){
		unsigned all_indices[ndims_spatial+1];
		for(unsigned i=0; i<ndims_spatial; i++) all_indices[i] = dir_ind[i];
		all_indices[ndims_spatial] = 0;
		unsigned base_ind = data.direct_index(all_indices);
		unsigned nbins = data.axes[ndims_spatial].size();
		for(unsigned i=0; i<nbins; i++){
			data.y0[base_ind+i] *= r;
		}
	}

	//--------------------------------------------------------------
	// MPI scatter the spectrum contents.
	// Must rescale zone stop list to account for number of groups
	//--------------------------------------------------------------
	void mpi_sum_scatter(vector<unsigned>& zone_stop_list){
		unsigned ngroups = data.axes[nuGridIndex].size();
		vector<unsigned> stop_list = zone_stop_list;
		for(unsigned i=0; i<stop_list.size(); i++) stop_list[i] *= ngroups;
		data.mpi_sum_scatter(stop_list);
	}
	void mpi_sum(){
		data.mpi_sum();
	}


	//--------------------------------------------------------------
	// Write data to specified location in an HDF5 file
	//--------------------------------------------------------------
	void write_hdf5_data(H5::H5File file, const string name) {
		data.write_HDF5(file, name);
	}

	//--------------------------------------------------------------
	// Write distribution function coordinates to an HDF5 file
	//--------------------------------------------------------------
	void write_hdf5_coordinates(H5::H5File /*file*/, const string /*name*/) const {
		// no extra axes for moment array
	}

	void add_isotropic(const unsigned dir_ind[NDIMS+1], const double E){
		unsigned indices[data.Ndims()];
		for(unsigned i=0; i<ndims_spatial; i++) indices[i] = dir_ind[i];
		indices[nuGridIndex] = dir_ind[ndims_spatial];

		// increment moments
		Tuple<double,4> tmp;
		tmp[0] = E;
		tmp[1] = 0;
		tmp[2] = E/3.;
		tmp[3] = 0;
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
