#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "grid_1D_sphere.h"
#include "physical_constants.h"
#include <fstream>
#include <iostream>

namespace pc = physical_constants;
using namespace std;


//------------------------------------------------------------
// initialize the zone geometry
//------------------------------------------------------------
void grid_1D_sphere::read_model_file(Lua* lua)
{
  // open up the model file, complaining if it fails to open
  string model_file = lua->scalar<string>("model_file");
  std::ifstream infile;
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

  // number of zones
  int n_zones;
  infile >> n_zones;
  z.resize(n_zones);
  r_out.resize(n_zones);
  vol.resize(n_zones);

  // other properties
  double texp;
  infile >> r_inner;
  infile >> texp;

  // read zone properties
  double r0;
  for(int i=0; i<n_zones; i++)
  {
    infile >> r_out[i];
    infile >> z[i].rho;
    infile >> z[i].T_gas;
    infile >> z[i].ni56;

    z[i].v[0] = r_out[i]/texp;
    z[i].e_rad = pc::a*pow(z[i].T_gas,4);
    if(i==0) r0 = r_inner;
    else r0 = r_out[i-1];
    vol[i] = 4.0*pc::pi/3.0*(r_out[i]*r_out[i]*r_out[i] - r0*r0*r0);
  }
  
  infile.close();
}

void grid_1D_sphere::custom_model(Lua* lua)
{

}

//------------------------------------------------------------
// Overly simple search to find zone
//------------------------------------------------------------
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


//------------------------------------------------------------
// return volume of zone (precomputed)
//------------------------------------------------------------
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

