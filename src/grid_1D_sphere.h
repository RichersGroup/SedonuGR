#ifndef _GRID_1D_SPHERE_H
#define _GRID_1D_SPHERE_H 1

#include "grid_general.h"
#include "locate_array.h"
#include <vector>

//*******************************************
// 1-Dimensional Spherical geometry
//*******************************************
class grid_1D_sphere: public grid_general
{

private:

  // store location of the outer edge of the zone.
  locate_array r_out;

  // store volumes explicitly
  std::vector<double> vol;

public:

  virtual ~grid_1D_sphere() {}

  void read_model_file(Lua* lua);
  void custom_model(Lua* lua);

  // required functions
  int    get_zone(double *);
  double zone_volume(int);
  double zone_min_length(int);
  void   sample_in_zone(int, std::vector<double>, double[3]);
  void   velocity_vector(int i, double[3], double[3]);
  void   write_ray(int iw);
  void   coordinates(int i,double r[3]) {
    r[0] = r_out[i]; r[1] = 0; r[2] = 0;}
};


#endif
