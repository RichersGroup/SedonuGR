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

//------------------------------------------------------------------
//*****************************************************************
//*************************  GRID ********************************
//*****************************************************************
// The grid class is a construct whose main purpose is to handle 
// the geometry of the model.  It does two main things: (1) Reads
// in the input density,temperature,composition files (that of 
// course must have a specific geometry). (2) Given a set of 
// 3-d coordinates, it will give the corosponding zone index (or
// note that the coords are off the grid).
//
// The grid holds an array of zones, where key fluid data is stored
//
// The grid class is an abstract class that will be used to
// create subclasses (e.g. grid_1D_sphere, grid_3D_cart, etc...
//*****************************************************************


#ifndef _GRID_GENERAL_H
#define _GRID_GENERAL_H 1
#include "global_options.h"
#include <iostream>
#include "zone.h"
#include "Lua.h"
#include "particle.h"
#include "transport.h"

class transport;

using namespace std;

class grid_general
{

protected:

	// fill the grid with data from a model file
	virtual void read_model_file(Lua* lua) = 0;

	// output options
	int output_distribution;
	int output_hdf5;
	int do_annihilation;

public:

	virtual ~grid_general() {}

	string grid_type;

	// vector of zones
	std::vector<zone> z;

#ifdef __INTEL_COMPILER
	static const double tiny = 1e-5; // used to overshoot boundary to account for error in boundary distance calculation
#else
	static constexpr double tiny = 1e-5; // used to overshoot boundary to account for error in boundary distance calculation
#endif

	// set everything up
	void init(Lua* lua);

	// write out zone information
	void write_zones(const int iw) const;
	void write_header(ofstream& outf) const;
	void write_line(ofstream& outf, const int z_ind) const;
	virtual void write_rays(const int iw) const = 0;
	void write_hdf5_data(H5::H5File file) const;
	virtual void write_hdf5_coordinates(H5::H5File file) const = 0;

	//****** virtual functions (geometry specific)

	// get directional indices from the zone index
	virtual void zone_directional_indices(const int z_ind, vector<int>& dir_ind) const = 0;
	virtual void dims(vector<hsize_t>& dims) const = 0;
	virtual hsize_t dimensionality() const = 0;

	// get the velocity squared from the stored velocity vector
	double zone_speed2(const int z_ind) const;

	// get zone index from x,y,z position
	virtual int zone_index(const vector<double>& x) const   = 0;

	// return volume of zone z_ind
	virtual double zone_lab_volume(const int z_ind) const         = 0;

	// return rest mass in cell
	double zone_rest_mass(const int z_ind) const;
	double zone_comoving_volume(const int z_ind) const;
	bool good_zone(const int z_ind) const;
	double total_rest_mass() const;

	// return the smallest length dimension of zone  z_ind
	virtual double zone_min_length(const int z_ind) const     = 0;

	// randomly sample a position within the zone z_ind
	virtual void cartesian_sample_in_zone(const int z_ind,const vector<double>& rand, vector<double>& x) const = 0;

	// give the velocity vector at this point
	virtual void cartesian_velocity_vector(const vector<double>& x, vector<double>& v, int z_ind=-1) const = 0;

	// get the coordinates at the center of the zone z_ind
	virtual void zone_coordinates(const int z_ind, vector<double>& r) const = 0;
	virtual double zone_radius(const int z_ind) const = 0;

	// boundary conditions
	virtual void reflect_outer(particle *p) const = 0;
	virtual double lab_dist_to_boundary(const particle *p) const = 0;
	virtual void symmetry_boundaries(particle* p) const = 0;
};


#endif

