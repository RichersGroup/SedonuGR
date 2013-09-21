#include <algorithm>
#include <mpi.h>
#include <math.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
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
  const int verbose = (my_rank == 0);

  // open up the model file, complaining if it fails to open
  string model_file = lua->scalar<string>("model_file");
  ifstream infile;
  infile.open(model_file.c_str());
  if(infile.fail())
  {
    cout << "Error: can't read the model file." << model_file << endl;
    exit(4);
  }

  // geometry of model
  infile >> grid_type;
  if(grid_type != "1D_sphere"){
    cout << "Error: grid_type parameter disagrees with the model file." << endl;
  }

  // type of system
  string system;
  infile >> system;

  // number of zones
  int n_zones;
  infile >> n_zones;
  z.resize(n_zones);
  r_out.resize(n_zones);
  vol.resize(n_zones);

  // other properties
  double texp;
  infile >> r_out.min;
  infile >> texp;

  // read zone properties for a supernova remnant
  double r0;
  if(system == "SNR") for(int i=0; i<n_zones; i++)
  {
    if(i==0 && verbose) cout << "# Reading SNR model file..." << endl;
    infile >> r_out[i];
    infile >> z[i].rho;
    infile >> z[i].T_gas;
    infile >> z[i].ni56;

    z[i].v[0] = r_out[i]/texp;
    z[i].e_rad = pc::a*pow(z[i].T_gas,4);
    if(i==0) r0 = r_out.min;
    else r0 = r_out[i-1];
    vol[i] = 4.0*pc::pi/3.0*(r_out[i]*r_out[i]*r_out[i] - r0*r0*r0);

    z[i].Ye = -1e99; // poison
  }

  // read zone properties for a gamma-ray burst
  else if(system == "GRB") for(int i=0; i<n_zones; i++)
  {
    if(i==0 && verbose) cout << "# Reading GRB model file..." << endl;
    infile >> r_out[i];
    infile >> z[i].rho;
    infile >> z[i].T_gas;
    infile >> z[i].Ye;

    z[i].v[0] = 0;
    z[i].e_rad = pc::a*pow(z[i].T_gas,4);
    if(i==0) r0 = r_out.min;
    else r0 = r_out[i-1];
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

//------------------------------------------------------------
// Overly simple search to find zone
//------------------------------------------------------------
int grid_1D_sphere::get_zone(double *x)
{
  double r = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);

  // check if off the boundaries
  if(r <= r_out.min             ) return -1;
  if(r >  r_out[r_out.size()-1] ) return -2;

  // find in zone array using stl algorithm upper_bound and subtracting iterators
  return r_out.locate(r);
}


//------------------------------------------------------------
// return volume of zone (precomputed)
//------------------------------------------------------------
double  grid_1D_sphere::zone_volume(int i)
{
  return vol[i];
}


//------------------------------------------------------------
// return length of zone
//------------------------------------------------------------
double  grid_1D_sphere::zone_min_length(int i)
{
  if (i == 0) return (r_out[i] - r_out.min);
  else return (r_out[i] - r_out[i-1]);
}



//------------------------------------------------------------
// sample a random position within the spherical shell
//------------------------------------------------------------
void grid_1D_sphere::sample_in_zone
(int i, std::vector<double> ran, double r[3])
{
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
void grid_1D_sphere::velocity_vector(int i, double x[3], double v[3])
{
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

