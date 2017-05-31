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

#include "LocateArray.h"
#include "Grid.h"
#include "H5Cpp.h"
#include <fstream>

class Grid;

class SpectrumArray {

public:

	virtual ~SpectrumArray() {}

	// MPI functions
	virtual void MPI_average(const int proc) = 0;

	// Count a packets
	virtual void count(const double D[3], const int Dsize, const double nu, const double E) = 0;

	//  void normalize();
	virtual void rescale(const double) = 0;
	virtual void wipe() = 0;

	// integrate over nu,mu,phi
	virtual double average_nu() const = 0;
	virtual double integrate() const = 0;
	virtual void integrate_over_direction(std::vector<double>& edens) const = 0;

	// Print out
	virtual void print(const int iw, const int species) const = 0;
	virtual void write_hdf5_data(H5::H5File file, const int s, const int dir_ind[], const hsize_t n_spatial_dims) const = 0;
	virtual void write_hdf5_coordinates(H5::H5File file, const Grid* grid) const = 0;
	virtual void write_header(std::ofstream& outf) const = 0;
	virtual void write_line(std::ofstream& outf) const = 0;

};

#endif
