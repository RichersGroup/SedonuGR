#ifndef _GRID_2D_SPHERE_H
#define _GRID_2D_SPHERE_H 1

#include "grid_general.h"
#include "locate_array.h"
#include <vector>
#include "global_options.h"

//*******************************************
// 1-Dimensional Spherical geometry
//*******************************************
class grid_2D_sphere: public grid_general
{

private:

	static const int dimensionality = 2;

	// store location of the outer edges of the zone.
	// order of zone array: r is increased fastest
	locate_array r_out;
	locate_array theta_out;

public:

	virtual ~grid_2D_sphere() {}

	void read_model_file(Lua* lua);
	void custom_model(Lua* lua);

	// required functions
	int    zone_index             (const vector<double>& x                                       ) const;
	int    zone_index             (const int i, const int j                                      ) const;
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
