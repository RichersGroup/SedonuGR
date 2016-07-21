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

#ifndef _GRID_2D_SPHERE_H
#define _GRID_2D_SPHERE_H 1

#include "GridNewtonian.h"

// PARAMETERS
//    Grid2DSphere_Flash_rgrid_file
//    Grid2DSphere_Flash_thetagrid_file
//    Grid2DSphere_Nagakura_xcoords_file
//    Grid2DSphere_Nagakura_ycoords_file

//*******************************************
// 1-Dimensional Spherical geometry
//*******************************************
class Grid2DSphere: public GridNewtonian
{

private:
	// store location of the outer edges of the zone.
	// order of zone array: r is increased fastest
	LocateArray r_out;
	LocateArray theta_out;

public:

	virtual ~Grid2DSphere() {}

	void read_model_file(Lua* lua);
	void read_flash_file(Lua* lua);
	void read_nagakura_file(Lua* lua);
	void custom_model(Lua* lua);

	// required functions
	int    zone_index             (const double x[3], const int xsize                            ) const;
	int    zone_index             (const int i, const int j                                      ) const;
	double zone_lab_volume        (const int z_ind                                               ) const;
	double zone_min_length        (const int z_ind                                               ) const;
	void zone_coordinates         (const int z_ind, double r[2], const int rsize                 ) const;
	void zone_directional_indices (const int z_ind, int dir_ind[2], const int size               ) const;
	void cartesian_sample_in_zone (const int z_ind, const double rand[3], const int randsize, double x[3], const int xsize) const;
	void cartesian_velocity_vector(const double x[3], const int xsize, double v[3], const int vsize, int z_ind) const;
	void write_rays               (const int iw                                                  ) const;
	void reflect_outer            (LorentzHelper *lh                                             ) const;
	void symmetry_boundaries      (LorentzHelper *lh                                             ) const;
	double lab_dist_to_boundary   (const LorentzHelper *lh                                       ) const;
	double zone_radius            (const int z_ind) const;
	void dims                     (hsize_t dims[2], const int size) const;
	hsize_t dimensionality() const {return 2;};
	void write_hdf5_coordinates(H5::H5File file) const;
};


#endif
