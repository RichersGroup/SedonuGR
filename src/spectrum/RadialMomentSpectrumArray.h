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
	MultiDArray<double,ndims_spatial+2> data;

	unsigned nuGridIndex, momGridIndex, nranks;

public:

	//--------------------------------------------------------------
	// Initialization and Allocation
	//--------------------------------------------------------------
	void init(const vector<Axis>& spatial_axes, const Axis& nu_grid, const int max_rank) {
		vector<Axis> axes;
		for(int i=0; i<spatial_axes.size(); i++) axes.push_back(spatial_axes[i]);

		axes.push_back(nu_grid);
		nuGridIndex = axes.size()-1;

		// moment axis
		nranks = max_rank+1;
		double minval = -0.5; // this is hacked...don't know a better way to include them all together
		vector<double> top(nranks);
		vector<double> mid(nranks);
		for(int i=0; i<nranks; i++){
			mid[i] = i;
			top[i] = i + 0.5;
		}
		axes.push_back(Axis(minval, top, mid));
		momGridIndex = axes.size()-1;

		// set up the data structure
		data = MultiDArray<double,ndims_spatial+2>(axes);
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
	void count(const double D[3], const vector<unsigned>& dir_ind, const double nu, const double E) {
		PRINT_ASSERT(E, >=, 0);
		PRINT_ASSERT(E, <, INFINITY);

		unsigned indices[data.Ndims()];
		for(int i=0; i<dir_ind.size(); i++) indices[i] = dir_ind[i];

		int nu_bin = data.axes[nuGridIndex].bin(nu);
		nu_bin = max(nu_bin, 0);
		nu_bin = min(nu_bin, data.axes[nuGridIndex].size()-1);
		indices[nuGridIndex] = nu_bin;

		// increment moments
		double tmp = E;
		for (int rank = 0; rank<nranks; rank++) {
			indices[momGridIndex] = rank;
			data.add(indices, tmp);
			tmp *= D[2];
		}
	}

	void rescale(double r) {
		for(unsigned i=0;i<data.size();i++) data.y0[i] *= r;
	}


	//--------------------------------------------------------------
	// MPI average the spectrum contents
	//--------------------------------------------------------------
	// only process 0 gets the reduced spectrum to print
	void MPI_average(){
		data.MPI_combine();
	}


	//--------------------------------------------------------------
	// Write data to specified location in an HDF5 file
	//--------------------------------------------------------------
	void write_hdf5_data(H5::H5File file, const string name) const {
		data.write_HDF5(file, name);
	}

	//--------------------------------------------------------------
	// Write distribution function coordinates to an HDF5 file
	//--------------------------------------------------------------
	void write_hdf5_coordinates(H5::H5File file, const string name) const {
		// no extra axes for moment array
	}
};

#endif
