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

#ifndef _GRID_0D_ISOTROPIC_H
#define _GRID_0D_ISOTROPIC_H 1

#include "Grid.h"

//*******************************************
// 0-Dimensional Isotropic geometry
//*******************************************
class Grid0DIsotropic: public Grid
{

public:

	Grid0DIsotropic();
	virtual ~Grid0DIsotropic() {}

	void read_model_file(Lua* lua);

	// required functions
	int  zone_index                (const Tuple<double,4>& x                                 ) const;
	double zone_lab_3volume        (const int z_ind                                   ) const;
	double zone_min_length         (const int z_ind                                   ) const;
	double d_boundary(const EinsteinHelper *eh) const;
	double d_randomwalk(const EinsteinHelper *eh) const;
	void zone_coordinates          (const int z_ind, double r[0], const int rsize     ) const;
	void zone_directional_indices  (const int z_ind, vector<unsigned>& dir_ind        ) const;
	Tuple<double,4> sample_in_zone (const int z_ind, ThreadRNG* rangen                ) const;
	void interpolate_fluid_velocity(EinsteinHelper* eh                                ) const;
	void symmetry_boundaries       (EinsteinHelper *eh                                ) const;
	double zone_radius             (const int z_ind                                   ) const;
	void dims                      (hsize_t dims[0], const int size                   ) const;
	double zone_lorentz_factor     (const int z_ind                                   ) const;
	hsize_t dimensionality() const {return 0;};
	void axis_vector(vector<Axis>& axes) const;
	void write_child_zones(H5::H5File file);

	// GR functions
	Tuple<double,4> dk_dlambda(const EinsteinHelper& eh) const; // Gamma^alhpa_mu_nu
	void interpolate_shift(EinsteinHelper* eh) const;
	void interpolate_3metric(EinsteinHelper* eh) const;
	void grid_coordinates(const Tuple<double,4>& xup, double coords[NDIMS]) const;
};


#endif
