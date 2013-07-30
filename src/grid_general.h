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

using namespace std;

class grid_general
{

 public:

  string grid_type;

  // vector of zones
  std::vector<zone> z;
  int n_zones;

  // mpi reduce quantities
  void reduce_radiation();
  void reduce_radiation_block(int, int);
  void reduce_T_gas();

  //****** virtual functions (geometry specific)

  // get zone index from x,y,z position
  virtual int get_zone(double *)   = 0;

  // return volume of zone i
  virtual double zone_volume(int i)         = 0;

  // return the smallest length dimension of zone  i
  virtual double zone_min_length(int i)     = 0;
  
  // randomly sample a position within the zone i
  virtual void sample_in_zone(int,std::vector<double>,double[3]) = 0;
  
  // give the velocity vector at this point in zone i
  virtual void velocity_vector(int i, double[3], double[3]) = 0;
  
  // get the coordinates at the center of the zone i
  virtual void coordinates(int i,double r[3]) = 0;

};


#endif

