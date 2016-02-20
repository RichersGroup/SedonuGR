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

#include "locate_array.h"
#include "H5Cpp.h"

// default values
#define DEFAULT_NAME "spectrum_array"

class spectrum_array {

private:

	// bin arrays
	// values represent bin upper walls (the single locate_array.min value is the leftmost wall)
	// underflow is combined into leftmost bin (right of the locate_array.min)
	// overflow is combined into the rightmost bin (left of locate_array[size-1])
	locate_array nu_grid;
	locate_array mu_grid;
	locate_array phi_grid;

	// counting arrays
	std::vector<double> flux;

public:

	// Initialize
	void init(const locate_array nu_grid, const locate_array mu_grid, const locate_array phi_grid);
	void init(const std::vector<double> nu_grid, const int n_mu, const int n_phi);

	// MPI functions
	void MPI_average();

	// Count a packets
	void count(const double D[3], const int Dsize, const double nu, const double E);

	//  void normalize();
	void rescale(const double);
	void wipe();

	// integrate over nu,mu,phi
	double average_nu() const;
	double integrate() const;
	void integrate_over_direction(std::vector<double>& edens) const;

	// get bin centers and indices
	unsigned size() const;
	double get(const int index) const;
	double  nu_center(const unsigned index) const;
	double  mu_center(const unsigned index) const;
	double phi_center(const unsigned index) const;
	double  nu_bin_center(const unsigned index) const;
	double  mu_bin_center(const unsigned index) const;
	double phi_bin_center(const unsigned index) const;

	// Print out
	void print(const int iw, const int species) const;
	void write_hdf5_data(H5::DataSet dataset, H5::DataSpace dataspace) const;
	void write_hdf5_coordinates(H5::H5File file) const;

	// Indexing
	unsigned index(const unsigned nu_bin,const unsigned mu_bin,const unsigned phi_bin) const;
	unsigned  nu_bin(const unsigned index) const;
	unsigned  mu_bin(const unsigned index) const;
	unsigned phi_bin(const unsigned index) const;
	unsigned  nu_dim() const {return  nu_grid.size();};
	unsigned  mu_dim() const {return  mu_grid.size();};
	unsigned phi_dim() const {return phi_grid.size();};
};

#endif
