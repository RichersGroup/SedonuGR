#ifndef _GRID_3D_CART_H
#define _GRID_3D_CART_H 1

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
  int *ix,*iy,*iz;
  int reflect_x, reflect_y, reflect_z;

public:

  virtual ~grid_3D_cart() {}

  void read_model_file(Lua* lua);
  void custom_model(Lua* lua);

  // required functions
  double zone_speed2(const int z_ind) const;
  int       zone_index(const double *) const;
  double    zone_volume(const int) const;
  double    zone_min_length(const int) const;
  void      sample_in_zone(const int, const std::vector<double>, double[3]) const;
  void      velocity_vector(const int i, const double[3], double[3]) const;
  void      print() const;
  void      cartesian_coordinates(const int z_ind, vector<double>& r) const;
  void      write_rays(const int iw) const;
  void reflect_outer(particle *) const;
  double dist_to_boundary(const particle *) const;
};


#endif
