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

	static const unsigned nelements = 6;
	MultiDArray<nelements,2> data;
	unsigned nuGridIndex;

	//--------------------------------------------------------------
	// Initialization and Allocation
	//--------------------------------------------------------------
	void init(const vector<Axis>& spatial_axes, const Axis& nu_grid) {
		vector<Axis> axes;
		for(int i=0; i<spatial_axes.size(); i++) axes.push_back(spatial_axes[i]);

		axes.push_back(nu_grid);
		nuGridIndex = axes.size()-1;

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
	void count(const double D[3], const unsigned dir_ind[1+1], const double nu, const double E) {
		PRINT_ASSERT(E, >=, 0);
		PRINT_ASSERT(E, <, INFINITY);

		unsigned base_ind = data.direct_index(dir_ind);

		Tuple<double,6> tmp = data[base_ind];
		// increment moments. Take advantage of fact that moments are stored in order.
		tmp[0] =  E;      // E		indices[momGridIndex] = 0;

		tmp[0] =  E*D[2]; // F
		tmp[0] =  E*D[2]*D[2]; // P^rr
		tmp[0] =  E * (D[0]*D[0] + D[1]*D[1])*0.5; // average of P^tt and P^pp
		tmp[0] =  E * D[2]*D[2]*D[2]; // W^rrr
		tmp[0] =  E * D[2]*(D[0]*D[0] + D[1]*D[1])*0.5; // average of W^rtt and W^rpp
		data.add(dir_ind,tmp);
	}


	void rescale(double r) {
		for(unsigned i=0;i<data.size();i++) data.y0[i] *= r;
	}


	//--------------------------------------------------------------
	// MPI average the spectrum contents
	//--------------------------------------------------------------
	// only process 0 gets the reduced spectrum to print
	void MPI_average()
	{
		data.MPI_AllCombine();
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
