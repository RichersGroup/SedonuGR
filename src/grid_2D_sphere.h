#ifndef _GRID_1D_SPHERE_H
#define _GRID_1D_SPHERE_H 1

#include "grid_general.h"
#include "locate_array.h"
#include <vector>

//*******************************************
// 1-Dimensional Spherical geometry
//*******************************************
class grid_2D_sphere: public grid_general
{

private:

  // store location of the outer edges of the zone.
  // order of zone array: r is increased fastest
  locate_array r_out;
  locate_array theta_out;

  // store volumes explicitly
  std::vector<double> vol;

public:

  virtual ~grid_1D_sphere() {}

  void read_model_file(Lua* lua);
  void custom_model(Lua* lua);

  // required functions
  double zone_speed2(const int z_ind) const;
  int    get_zone(double *);
  double zone_volume(int);
  double zone_min_length(int);
  void   sample_in_zone(int, std::vector<double>, double[3]);
  void   velocity_vector(int i, double[3], double[3]);
  void   write_ray(int iw);
  void   coordinates(int i,double r[3]);
  void reflect_outer(particle *) const;
  double dist_to_boundary(const particle *) const;
};


#endif
