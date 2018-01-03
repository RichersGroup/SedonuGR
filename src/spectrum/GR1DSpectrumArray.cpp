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
#include "GR1DSpectrumArray.h"
#include "Transport.h"

using namespace std;
namespace pc = physical_constants;


//--------------------------------------------------------------
// Initialization and Allocation
//--------------------------------------------------------------
void GR1DSpectrumArray::init(const vector<Axis>& spatial_axes, const Axis& nu_grid) {
	vector<Axis> axes;
	for(int i=0; i<spatial_axes.size(); i++) axes.push_back(spatial_axes[i]);

	axes.push_back(nu_grid);
	nuGridIndex = axes.size()-1;

	// moment axis
	double minval = 0;
	vector<double> top(nelements);
	vector<double> mid(nelements);
	for(int i=0; i<nelements; i++){
		mid[i] = i + 0.5;
		top[i] = i + 1;
	}
	axes.push_back(Axis(minval, top, mid));
	momGridIndex = axes.size()-1;
	
	// set up the data structure
	data = MultiDArray<2>(axes);
	data.wipe();
}

//--------------------------------------------------------------
// Functional procedure: Wipe
//--------------------------------------------------------------
void GR1DSpectrumArray::wipe() {
	data.wipe();
}


//--------------------------------------------------------------
// count a particle
////--------------------------------------------------------------
void GR1DSpectrumArray::count(const double D[3], const vector<unsigned>& dir_ind, const double nu, const double E) {
	PRINT_ASSERT(E, >=, 0);
	PRINT_ASSERT(E, <, INFINITY);

	unsigned indices[data.Ndims()];
	for(int i=0; i<dir_ind.size(); i++) indices[i] = dir_ind[i];

	int nu_bin = data.axes[nuGridIndex].bin(nu);
	nu_bin = max(nu_bin, 0);
	nu_bin = min(nu_bin, data.axes[nuGridIndex].size()-1);
	indices[nuGridIndex] = nu_bin;

	indices[momGridIndex] = 0;
	unsigned base_ind = data.direct_index(indices);

	// increment moments. Take advantage of fact that moments are stored in order.
	data.direct_add(base_ind  , E     ); // E
	data.direct_add(base_ind+1, E*D[2]); // E
	data.direct_add(base_ind+2, E*D[2]*D[2]); // P^rr
	data.direct_add(base_ind+3, E * (D[0]*D[0] + D[1]*D[1])*0.5); // average of P^tt and P^pp
	data.direct_add(base_ind+4, E * D[2]*D[2]*D[2]); // W^rrr
	data.direct_add(base_ind+5, E * D[2]*(D[0]*D[0] + D[1]*D[1])*0.5); // average of W^rtt and W^rpp
}


void GR1DSpectrumArray::rescale(double r) {
	for(unsigned i=0;i<data.size();i++) data.y0[i] *= r;
}


//--------------------------------------------------------------
// MPI average the spectrum contents
//--------------------------------------------------------------
// only process 0 gets the reduced spectrum to print
void GR1DSpectrumArray::MPI_average()
{
	data.MPI_AllCombine();
}


//--------------------------------------------------------------
// Write data to specified location in an HDF5 file
//--------------------------------------------------------------
void GR1DSpectrumArray::write_hdf5_data(H5::H5File file, const string name) const {
	data.write_HDF5(file, name);
}

//--------------------------------------------------------------
// Write distribution function coordinates to an HDF5 file
//--------------------------------------------------------------
void GR1DSpectrumArray::write_hdf5_coordinates(H5::H5File file, const string name) const {
	// no extra axes for moment array
}

