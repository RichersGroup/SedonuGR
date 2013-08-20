#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "grid_3D_cart.h"
#include "physical_constants.h"
#include <fstream>
#include <iostream>
#include "Lua.h"

namespace pc = physical_constants;

//************************************************************
// initialize the zone geometry
//************************************************************
void grid_3D_cart::custom_model(Lua* lua)
//(std::vector<int> n0, std::vector<double> ds, std::vector<double>s0)
{
  vector<int>n0 = lua->vector<int>("n0");
  vector<double>ds = lua->vector<double>("ds");
  vector<double>s0 = lua->vector<double>("s0");

  grid_type = "3D_cart";

  nx = n0[0];
  ny = n0[1];
  nz = n0[2];

  dx = ds[0];
  dy = ds[1];
  dz = ds[2];

  x0 = s0[0];
  y0 = s0[1];
  z0 = s0[2];

  vol = dx*dy*dz;
  if      ((dx < dy)&&(dx < dz)) min_ds = dx;
  else if ((dy < dx)&&(dy < dz)) min_ds = dx;
  else min_ds = dz;

  // allocate zones
  int n_zones = nx*ny*nz;
  z.resize(n_zones);

  // get a reverse map to x,y,z indices
  ix = new int[n_zones];
  iy = new int[n_zones];
  iz = new int[n_zones];
  int ind = 0;
  for (int i=0;i<nx;i++)
    for (int j=0;j<ny;j++)
      for (int k=0;k<nz;k++)
      {
	ix[ind] = i;
	iy[ind] = j;
	iz[ind] = k;
	ind++;
      }
}

void grid_3D_cart::read_model_file(Lua* lua)
{
  cout << "Error: there is no model file reader for grid type grid_3D_cart." << endl;
  exit(12);
}

//************************************************************
// Overly simple search to find zone
//************************************************************
int grid_3D_cart::get_zone(double *x)
{
  int i = floor(x[0]/dx - x0);
  int j = floor(x[1]/dy - y0);
  int k = floor(x[2]/dz - z0);

  // check for off grid
  if ((i < 0)||(i > nx-1)) return -2;
  if ((j < 0)||(j > nx-1)) return -2;
  if ((k < 0)||(k > nx-1)) return -2;
  
  int ind =  i*ny*nz + j*ny + k; 
  return ind;
}


//************************************************************
// return volume of zone (precomputed)
//************************************************************
double grid_3D_cart::zone_volume(int i)
{
  return vol;
}


//************************************************************
// sample a random position within the spherical shell
//************************************************************
void grid_3D_cart::sample_in_zone
(int i, std::vector<double> ran,double r[3])
{
  r[0] = x0 + (ix[i] + ran[0])*dx;
  r[1] = y0 + (iy[i] + ran[1])*dy;
  r[2] = z0 + (iz[i] + ran[2])*dz;
}


//************************************************************
// return length of zone
//************************************************************
double  grid_3D_cart::zone_min_length(int i)
{
  return min_ds;
}



//************************************************************
// get the velocity vector 
//************************************************************
void grid_3D_cart::velocity_vector(int i, double x[3], double v[3])
{
  // assuming radial velocity (may want to interpolate here)
  v[0] = z[i].v[0];
  v[1] = z[i].v[1];
  v[2] = z[i].v[2];
}


void grid_3D_cart::coordinates(int i,double r[3])
{
  r[0] = x0 + (ix[i]+0.5)*dx;
  r[1] = y0 + (iy[i]+0.5)*dy;
  r[2] = z0 + (iz[i]+0.5)*dz;

}


