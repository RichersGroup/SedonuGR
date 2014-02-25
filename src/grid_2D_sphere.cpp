#include <algorithm>
#include <mpi.h>
#include <math.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include "grid_2D_sphere.h"
#include "physical_constants.h"

namespace pc = physical_constants;
using namespace std;


//------------------------------------------------------------
// initialize the zone geometry from model file
//------------------------------------------------------------
void grid_2D_sphere::read_model_file(Lua* lua)
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
  if(grid_type != "2D_sphere"){
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
void grid_2D_sphere::custom_model(Lua* lua)
{
  cout << "Error: there is no custom model programmed for grid_2D_sphere." << endl;
  exit(11);
}

//------------------------------------------------------------
// Overly simple search to find zone
//------------------------------------------------------------
int grid_2D_sphere::get_zone(double *x)
{
  double r  = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
  double theta = atan( x[2] / sqrt(x[0]*x[0] + x[1]*x[1]) );

  // check if off the boundaries
  if(r     <= r_out.min                    ) return -1;
  if(r     >  r_out[r_out.size()-1]        ) return -2;
  if(theta <= theta_out.min                ) return -2;
  if(theta >  theta_out[theta_out.size()-1]) return -2;

  // find in zone array using stl algorithm upper_bound and subtracting iterators
  int i_r     =     r_out.locate(r    );
  int i_theta = theta_out.locate(theta);
  return i_r + i_theta*r_out.size();
}


//------------------------------------------------------------
// return volume of zone (precomputed)
//------------------------------------------------------------
double  grid_2D_sphere::zone_volume(int i)
{
  return vol[i];
}


//------------------------------------------------------------
// return length of zone
//------------------------------------------------------------
double  grid_2D_sphere::zone_min_length(int i)
{
  double r_len, theta_len;
  int i_theta = i / r_out.size();
  int i_r     = i % r_out.size();
  double r_in = i_r == 0 ? r_out.min : r_out[i_r-1];
  double theta_in = i_theta == 0 ? theta_out.min : theta_out[i_theta-1];

  // the 'minimum lengts' are just approximate.
  r_len     =     r_out[i_r    ] -     r_in;
  theta_len = theta_out[i_theta] - theta_in;
  
  // if r_in is zero, there will be problems, but simulations would not have done this.
  if(r_in == 0) return r_len;
  else return min(r_len, r_in*theta_len);
}



//------------------------------------------------------------
// sample a random position within the spherical shell
//------------------------------------------------------------
void grid_2D_sphere::sample_in_zone(int i, std::vector<double> ran, double r[3])
{
  // radius and theta indices
  int i_theta = i / r_out.size();
  int i_r     = i % r_out.size();

  // inner and outer radius of shell
  double r1 = i_r == 0 ? r_out.min : r_out[i_r-1];
  double r2 = r_out[i_r];

  // inner and outer cos(theta)
  double mu1 = i_theta == 0 ? theta_out.min : theta_out[i_theta-1];
  double mu2 = theta_out[i_theta];

  // sample radial position in shell using a probability integral transform
  double radius = pow( ran[0]*(r2*r2*r2 - r1*r1*r1) + r1*r1*r1, 1./3.);

  // sample cos(theta) uniformily
  double mu = mu1 + (mu2-mu1)*ran[1];
  double sin_theta = sqrt(1-mu*mu);

  // sample phi uniformily
  double phi = 2.0*pc::pi*ran[2];

  // set the real 3-d coordinates
  r[0] = radius*sin_theta*cos(phi);
  r[1] = radius*sin_theta*sin(phi);
  r[2] = radius*mu;
}



//------------------------------------------------------------
// get the velocity vector 
//------------------------------------------------------------
void grid_2D_sphere::velocity_vector(int i, double x[3], double v[3])
{
  // r and theta indices
  int i_theta = i / r_out.size();
  int i_r     = i % r_out.size();

  // may want to interpolate, but this is just using a step function
  v[0] = z[i].v[0];
  v[1] = z[i].v[1];
  v[2] = z[i]*v[2];
}


//------------------------------------------------------------
// Write the grid information out to a file
//------------------------------------------------------------
void grid_2D_sphere::write_ray(int iw)
{

}

//------------------------------------------------------------
// Return the (outer) coordinates of the cell
//------------------------------------------------------------
void grid_2D_sphere::coordinates(int i,double r[3]) {
  int i_theta = i / r_out.size();
  int i_r     = i % r_out.size();
  r[0] = r_out[i_r]; r[1] = r_out[i_theta]; r[2] = 0;
}
