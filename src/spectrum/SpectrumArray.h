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

#ifndef _SPECTRUM_ARRAY_H
#define _SPECTRUM_ARRAY_H 1

#include "H5Cpp.h"
#include "global_options.h"
#include <fstream>
#include <vector>
#include "EinsteinHelper.h"

using namespace std;

class SpectrumArray {

protected:

public:

	virtual ~SpectrumArray() {}

	// record some data
	virtual void count(const Tuple<double,4>& kup_tet, const unsigned dir_ind[NDIMS+1], const double E) = 0;

	// MPI functions
	virtual void mpi_sum_scatter(vector<unsigned>& zone_stop_list) = 0;
	virtual void mpi_sum() = 0;

	// Count a packets
	virtual void add_isotropic(const unsigned dir_ind[NDIMS+1], const double E) = 0;

	//  void normalize();
	virtual void rescale(const double) = 0;
	virtual void rescale_spatial_point(const unsigned dir_ind[NDIMS], const double) = 0;
	virtual void wipe() = 0;
	virtual double total() const = 0;

	// Print out
	virtual void write_hdf5_data(H5::H5File file, const string name) = 0;
	virtual void write_hdf5_coordinates(H5::H5File file, const string name) const = 0;
	virtual void annihilation_rate(const unsigned[] /*dir_ind[NDIMS]*/, const SpectrumArray* /*in_dist*/,
			const vector< vector<vector<double> > >& /*phi*/, const unsigned /*weight*/, Tuple<double,4>& /*fourforce*/) const{
		cout << "annihilation_rate is not implemented for this spectrum type!" << endl;
		assert(0);
	}
};

#endif
