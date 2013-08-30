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
  int       get_zone(double *);
  double    zone_volume(int);
  double    zone_min_length(int);
  void      sample_in_zone(int, std::vector<double>, double[3]);
  void      velocity_vector(int i, double[3], double[3]);
  void      print();
  void      coordinates(int i,double r[3]);

};


#endif
