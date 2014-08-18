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
  int    zone_index(const double *) const;
  double zone_volume(const int) const;
  double zone_min_length(const int) const;
  void   sample_in_zone(const int, const std::vector<double>, double[3]) const;
  void   velocity_vector(const int i, const double[3], double[3]) const;
  void   write_rays(const int iw) const;
  void   coordinates(const int i,double r[3]) const{
    r[0] = r_out[i]; r[1] = 0; r[2] = 0;}
  void reflect_outer(particle *) const;
  double dist_to_boundary(const particle *) const;
};


#endif
