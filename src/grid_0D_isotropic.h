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

#include "global_options.h"
#include "grid_general.h"
#include "locate_array.h"
#include <vector>

//*******************************************
// 0-Dimensional Isotropic geometry
//*******************************************
class grid_0D_isotropic: public grid_general
{

private:

	static const int dimensionality = 0;

public:

	virtual ~grid_0D_isotropic() {}

	void read_model_file(Lua* lua);
	void custom_model(Lua* lua);

	// required functions
	int  zone_index               (const vector<double>& x                                     ) const;
	double zone_lab_volume            (const int z_ind                                           ) const;
	double zone_min_length        (const int z_ind                                           ) const;
	void zone_coordinates         (const int z_ind, vector<double>& r                          ) const;
	void zone_directional_indices (const int z_ind, vector<int>& dir_ind                     ) const;
	void cartesian_sample_in_zone (const int z_ind, const vector<double>& rand, vector<double>& x) const;
	void cartesian_velocity_vector(const vector<double>& x, vector<double>& v, int z_ind     ) const;
	void write_rays               (const int iw                                              ) const;
	void reflect_outer            (particle *p                                               ) const;
	double lab_dist_to_boundary       (const particle *p                                         ) const;
	double zone_radius            (const int z_ind) const;
	void dims(vector<int>& dims) const;
};


#endif
