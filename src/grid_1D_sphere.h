#ifndef _GRID_1D_SPHERE_H
#define _GRID_1D_SPHERE_H 1

#include "grid_general.h"
#include "locate_array.h"
#include "global_options.h"
#include <vector>

//*******************************************
// 1-Dimensional Spherical geometry
//*******************************************
class grid_1D_sphere: public grid_general
{

private:

	static const int dimensionality = 1;

	// store location of the outer edge of the zone.
	locate_array r_out;


public:

	virtual ~grid_1D_sphere() {}

	void read_model_file(Lua* lua);
	void custom_model(Lua* lua);

	// required functions
	int  zone_index               (const vector<double>& x                                       ) const;
	double zone_lab_volume            (const int z_ind                                               ) const;
	double zone_min_length        (const int z_ind                                               ) const;
	void zone_coordinates         (const int z_ind, vector<double>& r                            ) const;
	void zone_directional_indices (const int z_ind, vector<int>& dir_ind                         ) const;
	void cartesian_sample_in_zone (const int z_ind, const vector<double>& rand, vector<double>& x) const;
	void cartesian_velocity_vector(const vector<double>& x, vector<double>& v                    ) const;
	void write_rays               (const int iw                                                  ) const;
	void reflect_outer            (particle *p                                                   ) const;
	double lab_dist_to_boundary       (const particle *p                                             ) const;
};


#endif
