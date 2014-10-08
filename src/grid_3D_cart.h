#ifndef _GRID_3D_CART_H
#define _GRID_3D_CART_H 1

#include "global_options.h"
#include "grid_general.h"
#include <vector>
#include "Lua.h"

//*******************************************
// 1-Dimensional Spherical geometry
//*******************************************
class grid_3D_cart: public grid_general
{

private:

	static const int dimensionality = 3;

	// specifics to this geometry

	int    nx, ny, nz; // number of zones in each dimension
	double dx, dy, dz; // length of each zone in each dimension
	double x0, y0, z0; // leftmost points
	double vol;        // volume of each zone = dx*dy*dz
	double min_ds;
	int reflect_x, reflect_y, reflect_z;

public:

	grid_3D_cart();
	virtual ~grid_3D_cart() {}

	void read_model_file(Lua* lua);
	void custom_model(Lua* lua);

	// required functions
	int    zone_index               (const vector<double>& x                                       ) const;
	int    zone_index               (const int i, const int j, const int k                         ) const;
	double zone_lab_volume              (const int z_ind                                               ) const;
	double zone_min_length          (const int z_ind                                               ) const;
	void   zone_coordinates         (const int z_ind, vector<double>& r                            ) const;
	void   zone_directional_indices (const int z_ind, vector<int>& dir_ind                         ) const;
	void   cartesian_sample_in_zone (const int z_ind, const vector<double>& rand, vector<double>& x) const;
	void   cartesian_velocity_vector(const vector<double>& x, vector<double>& v, int z_ind      ) const;
	void   write_rays               (const int iw                                                  ) const;
	void   reflect_outer            (particle *p                                                   ) const;
	double lab_dist_to_boundary         (const particle *p                                             ) const;
	double zone_radius              (const int z_ind) const;
};


#endif
