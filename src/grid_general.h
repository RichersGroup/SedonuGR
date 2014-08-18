//------------------------------------------------------------------
//*****************************************************************
//*************************  GRID ********************************
//*****************************************************************
// The grid class is a construct whose main purpose is to handle 
// the geometry of the model.  It does two main things: (1) Reads
// in the input density,temperature,composition files (that of 
// course must have a specific geometry). (2) Given a set of 
// 3-d coordinates, it will give the corosponding zone index (or
// note that the coords are off the grid).
//
// The grid holds an array of zones, where key fluid data is stored
//
// The grid class is an abstract class that will be used to
// create subclasses (e.g. grid_1D_sphere, grid_3D_cart, etc...
//*****************************************************************


#ifndef _GRID_GENERAL_H
#define _GRID_GENERAL_H 1
#include <string>
#include <iostream>
#include "zone.h"
#include "Lua.h"
#include "particle.h"

using namespace std;

class grid_general
{

 protected:

  // fill the grid with data from a model file
  virtual void read_model_file(Lua* lua) = 0;

  // fill the grid with data hard coded here
  virtual void custom_model(Lua* lua) = 0;

 public:

  virtual ~grid_general() {}

  string grid_type;

  // vector of zones
  std::vector<zone> z;

  // mpi reduce quantities
  /* void reduce_radiation(); */
  /* void reduce_radiation_block(int, int); */
  /* void reduce_T(); */
  /* void reduce_Ye(); */

  static const double tiny = 1e-3; // used to overshoot boundary to account for error in boundary distance calculation

  // set everything up
  void init(Lua* lua);

  // write out zone information
  void write_zones(const int iw) const;
  virtual void write_rays(const int iw) const = 0;

  //****** virtual functions (geometry specific)

  // get zone index from x,y,z position
  virtual int zone_index(const double *) const   = 0;

  // return volume of zone i
  virtual double zone_volume(const int i) const         = 0;

  // return the smallest length dimension of zone  i
  virtual double zone_min_length(const int i) const     = 0;
  
  // randomly sample a position within the zone i
  virtual void sample_in_zone(const int,const std::vector<double>,double[3]) const = 0;
  
  // give the velocity vector at this point in zone i
  virtual void velocity_vector(const int i, const double[3],double[3]) const = 0;
  
  // get the coordinates at the center of the zone i
  virtual void coordinates(const int i,double r[3]) const = 0;

  // boundary conditions
  virtual void reflect_outer(particle *) const = 0;
  virtual double dist_to_boundary(const particle *) const = 0;
};


#endif

