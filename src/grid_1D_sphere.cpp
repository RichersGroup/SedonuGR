#include <algorithm>
#include <mpi.h>
#include <math.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <cassert>
#include <limits>
#include "grid_1D_sphere.h"
#include "physical_constants.h"

namespace pc = physical_constants;
using namespace std;


//------------------------------------------------------------
// initialize the zone geometry from model file
//------------------------------------------------------------
void grid_1D_sphere::read_model_file(Lua* lua)
{
  // verbocity
  int my_rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &my_rank );
  const int rank0 = (my_rank == 0);
  if(rank0) cout << "# Reading 1D model file" << endl;

  // open up the model file, complaining if it fails to open
  string model_file = lua->scalar<string>("model_file");
  ifstream infile;
  infile.open(model_file.c_str());
  if(infile.fail()){
    cout << "Error: can't read the model file." << model_file << endl;
    exit(4);
  }

  // geometry of model
  infile >> grid_type;
  if(grid_type != "1D_sphere"){
    cout << "Error: grid_type parameter disagrees with the model file." << endl;
  }

  // number of zones
  int n_zones;
  infile >> n_zones;
  z.resize(n_zones,zone(1));
  r_out.resize(n_zones);
  vol.resize(n_zones);

  // read zone properties
  infile >> r_out.min;
  for(int i=0; i<n_zones; i++)
  {
    infile >> r_out[i];
    infile >> z[i].rho;
    infile >> z[i].T_gas;
    infile >> z[i].Ye;

    z[i].v[0] = 0;
    double r0 = r_out.bottom(i);
    vol[i] = 4.0*pc::pi/3.0*(r_out[i]*r_out[i]*r_out[i] - r0*r0*r0);
  }

  infile.close();
}

//------------------------------------------------------------
// Write a custom model here if you like
//------------------------------------------------------------
void grid_1D_sphere::custom_model(Lua* lua)
{
  cout << "Error: there is no custom model programmed for grid_1D_sphere." << endl;
  exit(11);
}

//------------------------------------
// get the velocity squared of a zone
//------------------------------------
double grid_1D_sphere::zone_speed2(const int z_ind) const{
	return z[z_ind].v[0];
}

//------------------------------------------------------------
// Overly simple search to find zone
//------------------------------------------------------------
int grid_1D_sphere::zone_index(const double *x) const
{
  double r = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);

  // check if off the boundaries
  if(r < r_out.min             ) return -1;
  if(r > r_out[r_out.size()-1] ) return -2;

  // find in zone array using stl algorithm upper_bound and subtracting iterators
  return r_out.locate(r);
}


//------------------------------------------------------------
// return volume of zone (precomputed)
//------------------------------------------------------------
double  grid_1D_sphere::zone_volume(const int i) const
{
	assert(i >= 0);
	return vol[i];
}


//------------------------------------------------------------
// return length of zone
//------------------------------------------------------------
double  grid_1D_sphere::zone_min_length(const int i) const
{
	assert(i >= 0);
	if (i == 0) return (r_out[i] - r_out.min);
	else return (r_out[i] - r_out[i-1]);
}



//------------------------------------------------------------
// sample a random position within the spherical shell
//------------------------------------------------------------
void grid_1D_sphere::sample_in_zone
(const int i, const std::vector<double> ran, double r[3]) const
{
	assert(i >= 0);
  // inner radius of shell
  double r1;
  if (i == 0) r1 = r_out.min;
  else r1 = r_out[i-1];

  // outer radius of shell
  double r2 = r_out[i];

  // sample radial position in shell using a probability integral transform
  double radius = pow( ran[0]*(r2*r2*r2 - r1*r1*r1) + r1*r1*r1, 1./3.);

  // random spatial angles
  double mu  = 1 - 2.0*ran[1];
  double phi = 2.0*pc::pi*ran[2];
  double sin_theta = sqrt(1 - mu*mu);

  // set the real 3-d coordinates
  r[0] = radius*sin_theta*cos(phi);
  r[1] = radius*sin_theta*sin(phi);
  r[2] = radius*mu;
}



//------------------------------------------------------------
// get the velocity vector 
//------------------------------------------------------------
void grid_1D_sphere::velocity_vector(const int i, const double x[3], double v[3]) const
{
	assert(i >= 0);
  // radius in zone
  double r = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);

  // assuming radial velocity (may want to interpolate here)
  // (the other two components are ignored and mean nothing)
  v[0] = x[0]/r*z[i].v[0];
  v[1] = x[1]/r*z[i].v[0];
  v[2] = x[2]/r*z[i].v[0];

  // check for pathological case
  if (r == 0)
  {
    v[0] = 0;
    v[1] = 0;
    v[2] = 0;
  }
}


//------------------------------------------------------------
// Write the grid information out to a file
//------------------------------------------------------------
void grid_1D_sphere::write_rays(const int iw) const
{
  // this is a 1D grid, so the function is exactly the same
  // as write_zones
}


//------------------------------------------------------------
// Reflect off the outer boundary
//------------------------------------------------------------
void grid_1D_sphere::reflect_outer(particle *p) const{
  double dr = 0;
  if(r_out.size()>1) dr = r_out[r_out.size()-1]-r_out[r_out.size()-2];
  else dr = r_out[0]-r_out.min;
  assert( fabs(p->r() - r_out[r_out.size()-1]) < tiny*dr);
  double velDotRhat = p->mu();
  double R = p->r();

  // invert the radial component of the velocity
  p->D[0] -= 2.*velDotRhat * p->x[0]/R;
  p->D[1] -= 2.*velDotRhat * p->x[1]/R;
  p->D[2] -= 2.*velDotRhat * p->x[2]/R;

  // put the particle just inside the boundary
  double newR = r_out[r_out.size()-1] - tiny*dr;
  p->x[0] = p->x[0]/R*newR;
  p->x[1] = p->x[1]/R*newR;
  p->x[2] = p->x[2]/R*newR;
  
  // must be inside the boundary, or will get flagged as escaped
  p->ind = zone_index(p->x);
  assert(p->r() < r_out[r_out.size()-1]);
}

//------------------------------------------------------------
// Find distance to outer boundary (less a tiny bit)
// negative distance means inner boundary
//------------------------------------------------------------
double grid_1D_sphere::dist_to_boundary(const particle *p) const{
  // Theta = angle between radius vector and direction (Pi if outgoing)
  // Phi   = Pi - Theta (angle on the triangle) (0 if outgoing)
  double Rout  = r_out[r_out.size()-1];
  double Rin   = r_out.min;
  double r  = p->r();
  double mu = p->mu();
  double d_outer_boundary = numeric_limits<double>::infinity();
  double d_inner_boundary = numeric_limits<double>::infinity();
  assert(r<Rout);
  assert(p->ind >= -1);

  // distance to inner boundary
  if(r >= Rin){
	  double radical = r*r*(mu*mu-1.0) + Rin*Rin;
	  if(Rin>0 && mu<0 && radical>=0){
		  d_inner_boundary = -r*mu - sqrt(radical);
		  assert(d_inner_boundary <= sqrt(Rout*Rout-Rin*Rin)*(1.0+tiny));
	  }
  }
  else{
	  d_inner_boundary = -r*mu + sqrt(r*r*(mu*mu-1.0) + Rin*Rin);
	  assert(d_inner_boundary <= 2.*Rin);
  }
  assert(d_inner_boundary >= 0);

  // distance to outer boundary
  d_outer_boundary = -r*mu + sqrt(r*r*(mu*mu-1.0) + Rout*Rout);
  assert(d_outer_boundary >= 0);
  assert(d_outer_boundary <= 2.*Rout);

  // make sure the particle ends up in a reasonable place
  return min(d_inner_boundary, d_outer_boundary);
}
