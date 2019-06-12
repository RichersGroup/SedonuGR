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

#ifndef _GR1D_SPECTRUM_ARRAY_H
#define _GR1D_SPECTRUM_ARRAY_H 1

#include "SpectrumArray.h"
#include "H5Cpp.h"
#include "global_options.h"
#include <fstream>
#include <mpi.h>
#include <sstream>
#include "Transport.h"

using namespace std;

class GR1DSpectrumArray : public SpectrumArray {

private:

	// bin arrays
	// values represent bin upper walls (the single locate_array.min value is the leftmost wall)
	// underflow is combined into leftmost bin (right of the locate_array.min)
	// overflow is combined into the rightmost bin (left of locate_array[size-1])

public:

	static const size_t nelements = 6;
	MultiDArray<double,nelements,2> data;

	//--------------------------------------------------------------
	// Initialization and Allocation
	//--------------------------------------------------------------
	void init(const vector<Axis>& spatial_axes, const Axis& nu_grid) {
		vector<Axis> axes;
		for(size_t i=0; i<spatial_axes.size(); i++) axes.push_back(spatial_axes[i]);

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
	void count(const Tuple<double,4>& kup_tet, const size_t dir_ind[NDIMS+1], const double E) {
		Tuple<double,3> D;
		for(size_t i=0; i<3; i++) D[i] = kup_tet[i];
		Metric::normalize_Minkowski<3>(D);
		PRINT_ASSERT(E, >=, 0);
		PRINT_ASSERT(E, <, INFINITY);

		Tuple<double,6> tmp;
		// increment moments. Take advantage of fact that moments are stored in order.
		tmp[0] =  E;      // E		indices[momGridIndex] = 0;
		tmp[1] =  E*D[2]; // F
		tmp[2] =  E*D[2]*D[2]; // P^rr
		tmp[3] =  E * (D[0]*D[0] + D[1]*D[1])*0.5; // average of P^tt and P^pp
		tmp[4] =  E * D[2]*D[2]*D[2]; // W^rrr
		tmp[5] =  E * D[2]*(D[0]*D[0] + D[1]*D[1])*0.5; // average of W^rtt and W^rpp
		data.add(dir_ind,tmp);
	}


	void rescale(double r) {
		for(size_t i=0;i<data.size();i++) data.y0[i] *= r;
	}
	void rescale_spatial_point(const size_t dir_ind[1], const double r){
		size_t all_indices[1+1];
		all_indices[0] = dir_ind[0];
		all_indices[1] = 0;
		size_t base_ind = data.direct_index(all_indices);
		size_t nbins = data.axes[1].size();
		for(size_t i=0; i<nbins; i++){
			data.y0[base_ind+i] *= r;
		}
	}


	//--------------------------------------------------------------
	// MPI scatter the spectrum contents.
	// Must rescale zone stop list to account for number of groups
	//--------------------------------------------------------------
	void mpi_sum_scatter(vector<size_t>& zone_stop_list){
		size_t ngroups = data.axes[1].size();
		vector<size_t> stop_list = zone_stop_list;
		for(size_t i=0; i<stop_list.size(); i++) stop_list[i] *= ngroups;
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
	void read_hdf5_data(H5::H5File file, const string name, const string /*axis_base*/) {
		vector<Axis> axes(2);
		axes[0].read_HDF5("/axes/x0(cm)",file);
		axes[1].read_HDF5("/axes/frequency(Hz)",file);

		data.read_HDF5(file, name, axes);
	}

	//--------------------------------------------------------------
	// Write distribution function coordinates to an HDF5 file
	//--------------------------------------------------------------
	void write_hdf5_coordinates(H5::H5File /*file*/, const string /*name*/) const {
		// no extra axes for moment array
	}

	void add_isotropic(const size_t dir_ind[NDIMS+1], const double E){
		PRINT_ASSERT(E, >=, 0);
		PRINT_ASSERT(E, <, INFINITY);

		Tuple<double,6> tmp;
		tmp = 0;

		// increment moments. Take advantage of fact that moments are stored in order.
		tmp[0] =  E;      // E		indices[momGridIndex] = 0;
		tmp[2] =  E/3.; // P^rr
		tmp[3] =  E/3.; // average of P^tt and P^pp
		data.add(dir_ind,tmp);
	}
	double total() const{
		double result=0;
		for(size_t i=0; i<data.size(); i++)
			result += data[i][0];
		return result;
	}
};

#endif
