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

#ifndef _GRID_3D_CART_H
#define _GRID_3D_CART_H 1

#include "GridNewtonian.h"

// PARAMETERS
//    Grid3DCart_THC_reflevel -- refinement level to use when reading in THC data
//    Grid3DCart_reflect_{x,y,z} -- reflect particles off the {x,y,z}=0 boundary. Also truncates the fluid grid there.
//    Grid3DCart_rotate_hemisphere_{x,y} -- assume 180-degree rotational symmetry in hemispheres where everywhere {x,y}>0
//    Grid3DCart_rotate_quadrant -- assume 90-degree rotational symmetry. x and y >0 everywhere.

//*******************************************
// 1-Dimensional Spherical geometry
//*******************************************
class Grid3DCart: public GridNewtonian
{

private:
	// specifics to this geometry

	int    nx[3]; // number of zones in each dimension
	double dx[3]; // length of each zone in each dimension
	double x0[3]; // leftmost points
	double x1[3]; // next-to leftmost points
	double xmax[3];
	int    reflect[3];
	int    rotate_hemisphere[2];
	int    rotate_quadrant;

public:

	Grid3DCart();
	virtual ~Grid3DCart() {}

	void get_deltas(const int z_ind, double delta[3], const int size) const;

	void read_model_file(Lua* lua);
	void read_SpEC_file(Lua* lua);
	void read_THC_file(Lua* lua);

	double zone_left_boundary(const unsigned dir, const unsigned dir_ind) const;
	double zone_right_boundary(const unsigned dir, const unsigned dir_ind) const;

	// required functions
	int    zone_index               (const double x[3], const int xsize                            ) const;
	int    zone_index               (const int i, const int j, const int k                         ) const;
	double zone_lab_volume          (const int z_ind                                               ) const;
	double zone_min_length          (const int z_ind                                               ) const;
	void   zone_coordinates         (const int z_ind, double r[3], const int rsize                 ) const;
	void   zone_directional_indices (const int z_ind, int dir_ind[3], const int size               ) const;
	void   sample_in_zone (const int z_ind, const double rand[3], const int randsize, double x[3], const int xsize) const;
	void   interpolate_fluid_velocity(const double x[3], const int sxize, double v[3], const int vsize, int z_ind) const;
	void   write_rays               (const int iw                                                  ) const;
	void   symmetry_boundaries      (LorentzHelper *lh                                             ) const;
	double lab_dist_to_boundary     (const LorentzHelper *lh                                       ) const;
	double zone_radius              (const int z_ind) const;
	void dims                       (hsize_t dims[3], const int size) const;
	hsize_t dimensionality() const {return 3;};
	void write_hdf5_coordinates(H5::H5File file) const;
	double zone_cell_dist(const double x_up[3], const int z_ind) const;
};


#endif
