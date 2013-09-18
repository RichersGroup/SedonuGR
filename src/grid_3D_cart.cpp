#include <math.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include "Lua.h"
#include "grid_3D_cart.h"
#include "physical_constants.h"

namespace pc = physical_constants;

//------------------------------------------------------------
// initialize the zone geometry
//------------------------------------------------------------
void grid_3D_cart::custom_model(Lua* lua)
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

//------------------------------------------------------------
// Read in a cartesian model file
//------------------------------------------------------------
void grid_3D_cart::read_model_file(Lua* lua)
{
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
  if(grid_type != "3D_cart"){
	cout << "Error: grid_type parameter disagrees with the model file." << endl;
	exit(4);
  }

  // type of system
  string system;
  infile >> system;

  // number of zones
  infile >> nx;
  infile >> ny;
  infile >> nz;
  // assume reflecting symmetry?
  infile >> reflect_x;
  infile >> reflect_y;
  infile >> reflect_z;
  if(reflect_x) nx*=2;
  if(reflect_y) ny*=2;
  if(reflect_z) nz*=2;
  int n_zones = nx*ny*nz;

  //set the zone sizes and volumes
  double xmax, ymax, zmax;
  infile >> x0;
  infile >> xmax;
  infile >> y0;
  infile >> ymax;
  infile >> z0;
  infile >> zmax;
  dx = (xmax-x0)/(double)nx; if(reflect_x) dx*=2;
  dy = (ymax-y0)/(double)ny; if(reflect_y) dy*=2;
  dz = (zmax-z0)/(double)nz; if(reflect_z) dz*=2;
  if      ((dx < dy)&&(dx < dz)) min_ds = dx;
  else if ((dy < dx)&&(dy < dz)) min_ds = dx;
  else min_ds = dz;
  vol = dx*dy*dz;

  // First loop - set indices and read zone values in file
  // loop order is the file order
  // index order matches that in get_zone
  ix = new int[n_zones];
  iy = new int[n_zones];
  iz = new int[n_zones];
  z.resize(n_zones);
  int ind = 0;
  bool rx,ry,rz;
  for (int k=0;k<nz;k++)
    for (int j=0;j<ny;j++)
      for (int i=0;i<nx;i++)
      {
    	// set current index
    	ind = i*ny*nz + j*nz + k;

    	// create reverse map to x,y,z indices
    	rx=false; ry=false; rz=false;
    	ix[ind] = i; if(reflect_x && i<nx/2) rx = true;
    	iy[ind] = j; if(reflect_y && j<ny/2) ry = true;
    	iz[ind] = k; if(reflect_z && k<nz/2) rz = true;

    	// read in values if not in reflected zone
    	if(!rx && !ry && !rz)
    	{
    	  infile >> z[ind].rho;
    	  infile >> z[ind].T_gas;
    	  infile >> z[ind].Ye;
    	  infile >> z[ind].v[0];
    	  infile >> z[ind].v[1];
    	  infile >> z[ind].v[2];
    	  // convert to proper units
    	  z[ind].T_gas /= (1e-6*pc::k_ev);
    	  z[ind].v[0] *= pc::c;
    	  z[ind].v[1] *= pc::c;
    	  z[ind].v[2] *= pc::c;
    	}
    	else{ //poison values
    	  double poison = -1e100;
      	  z[ind].rho   = poison;
      	  z[ind].T_gas = poison;
      	  z[ind].Ye    = poison;
      	  z[ind].v[0]  = poison;
      	  z[ind].v[1]  = poison;
      	  z[ind].v[2]  = poison;
    	}
      }

  // Second loop - apply symmetries
  ind=0;
  int origin_ind, origin_i, origin_j, origin_k;
  for (int k=0;k<nz;k++)
    for (int j=0;j<ny;j++)
      for (int i=0;i<nx;i++)
      {
    	  // set current index
    	  ind = i*ny*nz + j*nz + k;

    	  // are we in a reflection zone?
    	  rx=false; ry=false; rz=false;
    	  if(reflect_x && i<nx/2) rx = 1;
    	  if(reflect_y && j<ny/2) ry = 1;
    	  if(reflect_z && k<nz/2) rz = 1;

    	  // copy appropriate values
    	  if(rx) origin_i = (nx-1)-i; else origin_i = i;
    	  if(ry) origin_j = (ny-1)-j; else origin_j = j;
    	  if(rz) origin_k = (nz-1)-k; else origin_k = k;
    	  if(rx || ry || rz)
    	  {
    		  origin_ind = origin_i*ny*nz + origin_j*nz + origin_k;
    		  z[ind].rho   = z[origin_ind].rho;
    		  z[ind].T_gas = z[origin_ind].T_gas;
    		  z[ind].Ye    = z[origin_ind].Ye;
    		  z[ind].v[0]  = z[origin_ind].v[0];
    		  z[ind].v[1]  = z[origin_ind].v[1];
    		  z[ind].v[2]  = z[origin_ind].v[2];
    	  }
      }

  // adjust x0,y0,z0 to indicate new, reflected lower boundary
  cout << "zmin,zmax before adjusting:" << z0 << " " << zmax << endl;
  if(reflect_x) x0 = x0 - (xmax-x0);
  if(reflect_y) y0 = y0 - (ymax-y0);
  if(reflect_z) z0 = z0 - (zmax-z0);

  // do some output
  cout << "nx=" << nx << endl << "ny=" << ny << endl << "nz=" << nz << endl;
  cout << "number of zones:" << z.size() << endl;
  cout << "minima:{" << x0 << ", " << y0 << ", " << z0 << "}" << endl;
  cout << "maxima:{" << x0+(nx*dx) << ", " << y0+(ny*dy) << ", " << z0+(nz*dz) << "}" << endl;
  cout << "deltas:{" << dx << ", " << dy << ", " << dz << "}" << endl;
}

//------------------------------------------------------------
// Overly simple search to find zone
//------------------------------------------------------------
int grid_3D_cart::get_zone(double *x)
{
  int i = floor((x[0]-x0)/dx);
  int j = floor((x[1]-y0)/dy);
  int k = floor((x[2]-z0)/dz);

  // check for off grid
  if ((i < 0)||(i > nx-1)) return -2;
  if ((j < 0)||(j > ny-1)) return -2;
  if ((k < 0)||(k > nz-1)) return -2;
  
  int ind =  i*ny*nz + j*nz + k;
  return ind;
}


//------------------------------------------------------------
// return volume of zone (precomputed)
//------------------------------------------------------------
double grid_3D_cart::zone_volume(int i)
{
  return vol;
}


//------------------------------------------------------------
// sample a random position within the spherical shell
//------------------------------------------------------------
void grid_3D_cart::sample_in_zone
(int i, std::vector<double> ran,double r[3])
{
  r[0] = x0 + (ix[i] + ran[0])*dx;
  r[1] = y0 + (iy[i] + ran[1])*dy;
  r[2] = z0 + (iz[i] + ran[2])*dz;
}


//------------------------------------------------------------
// return length of zone
//------------------------------------------------------------
double  grid_3D_cart::zone_min_length(int i)
{
  return min_ds;
}



//------------------------------------------------------------
// get the velocity vector 
//------------------------------------------------------------
void grid_3D_cart::velocity_vector(int i, double x[3], double v[3])
{
  // may want to interpolate here
  v[0] = z[i].v[0];
  v[1] = z[i].v[1];
  v[2] = z[i].v[2];
}

//------------------------------------------------------------
// cell-centered coordinates of zone i
//------------------------------------------------------------
void grid_3D_cart::coordinates(int i,double r[3])
{
  r[0] = x0 + (ix[i]+0.5)*dx;
  r[1] = y0 + (iy[i]+0.5)*dy;
  r[2] = z0 + (iz[i]+0.5)*dz;
}
