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
#include <fstream>

using namespace std;

class GR1DSpectrumArray : public SpectrumArray {

private:

	// bin arrays
	// values represent bin upper walls (the single locate_array.min value is the leftmost wall)
	// underflow is combined into leftmost bin (right of the locate_array.min)
	// overflow is combined into the rightmost bin (left of locate_array[size-1])
	Axis nu_grid;

public:

	// counting arrays
	vector<double> moments;

	// Initialize
	void init(const Axis& nu_grid);

	// MPI functions
	virtual void MPI_average(const int proc);

	// Count a packets
	virtual void count(const double D[3], const double nu, const double E);

	//  void normalize();
	virtual void rescale(const double);
	virtual void wipe();

	// integrate over nu,mu,phi
	virtual double average_nu() const;
	virtual double integrate() const;
	virtual void integrate_over_direction(vector<double>& edens) const;

	// Print out
	virtual void print(const int iw, const int species) const;
	virtual void write_hdf5_data(H5::H5File file, const int s, const vector<unsigned>& dir_ind) const;
	virtual void write_hdf5_coordinates(H5::H5File file, const Grid* grid) const;
	virtual void write_header(ofstream& outf) const;
	virtual void write_line(ofstream& outf) const;

	// energy bin info
	double nu_bin_center(const unsigned index) const;

	// return sizes
	unsigned index(const int nu_bin, const int rank) const;
};

#endif
