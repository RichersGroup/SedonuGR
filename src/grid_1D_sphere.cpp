#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "grid_1D_sphere.h"
#include "physical_constants.h"
#include <fstream>
#include <iostream>

namespace pc = physical_constants;

//************************************************************
// initialize the zone geometry
//************************************************************
void grid_1D_sphere::init(double ri, std::vector<double> rr)
{
  grid_type = "1D_sphere";

  r_inner = ri;
  for (int i=0;i<rr.size();i++) 
  {
    r_out.push_back(rr[i]);

    double r0;
    if (i==0) r0 = r_inner; else r0 = rr[i-1];
    double v = 4.0*pc::pi/3.0*(rr[i]*rr[i]*rr[i] - r0*r0*r0);
    vol.push_back(v);
  }

  // allocate zones
  z.resize(r_out.size());
  n_zones = z.size();

}

//************************************************************
// Overly simple search to find zone
//************************************************************
int grid_1D_sphere::get_zone(double *x)
{
  double r = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);

  // check if inside inner edge
  if (r < r_inner) return -1;

  // find in zone array
  for (long int i=0;i<r_out.size();i++) 
    if (r < r_out[i]) return i;

  // off outer edge
  return -2;
}


//************************************************************
// return volume of zone (precomputed)
//************************************************************
double  grid_1D_sphere::zone_volume(int i)
{
  return vol[i];
}


//************************************************************
// return length of zone
//************************************************************
double  grid_1D_sphere::zone_min_length(int i)
{
  if (i == 0) return (r_out[i] - r_inner);
  else return (r_out[i] - r_out[i-1]);
}



//************************************************************
// sample a random position within the spherical shell
//************************************************************
void grid_1D_sphere::sample_in_zone
(int i, std::vector<double> ran, double r[3])
{
  // inner radius of shell
  double r_0;
  if (i == 0) r_0 = r_inner; 
  else r_0 = r_out[i-1];

  // thickness of shell
  double dr = r_out[i] - r_0;
  // sample radial position in shell
  r_0 = r_0 + dr*ran[0];

  // random spatial angles
  double mu  = 1 - 2.0*ran[1];
  double phi = 2.0*pc::pi*ran[2];
  double sin_theta = sqrt(1 - mu*mu);

  // set the real 3-d coordinates
  r[0] = r_0*sin_theta*cos(phi);
  r[1] = r_0*sin_theta*sin(phi);
  r[2] = r_0*mu;
}



//************************************************************
// get the velocity vector 
//************************************************************
void grid_1D_sphere::velocity_vector(int i, double x[3], double v[3])
{
  // radius in zone
  double rr = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);

  // assuming radial velocity (may want to interpolate here)
  v[0] = x[0]/rr*z[i].v[0];
  v[1] = x[1]/rr*z[i].v[0];
  v[2] = x[2]/rr*z[i].v[0];

  // check for pathological case
  if (rr == 0)
  {
    v[0] = 0;
    v[1] = 0;
    v[2] = 0;
  }
}

