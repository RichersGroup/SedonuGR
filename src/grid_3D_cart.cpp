#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cassert>
#include "mpi.h"
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
  // get mpi rank
  int my_rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &my_rank  );

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
  else if ((dy < dx)&&(dy < dz)) min_ds = dy;
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
  if(my_rank==0) cout << "# (zmin,zmax) before adjusting: (" << z0 << "," << zmax << ")" << endl;
  if(reflect_x) x0 = x0 - (xmax-x0);
  if(reflect_y) y0 = y0 - (ymax-y0);
  if(reflect_z) z0 = z0 - (zmax-z0);

  // debugging some output
  if(my_rank==0){
    cout << "# nx=" << nx << endl << "# ny=" << ny << endl << "# nz=" << nz << endl;
    cout << "# number of zones:" << z.size() << endl;
    cout << "# minima:{" << x0 << ", " << y0 << ", " << z0 << "}" << endl;
    cout << "# maxima:{" << x0+(nx*dx) << ", " << y0+(ny*dy) << ", " << z0+(nz*dz) << "}" << endl;
    cout << "# deltas:{" << dx << ", " << dy << ", " << dz << "}" << endl;
  }
}

//------------------------------------------------------------
// Overly simple search to find zone
//------------------------------------------------------------
int grid_3D_cart::get_zone(const double *x) const
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
double grid_3D_cart::zone_volume(const int i) const
{
  return vol;
}


//------------------------------------------------------------
// sample a random position within the cubical cell
//------------------------------------------------------------
void grid_3D_cart::sample_in_zone
(const int i, const std::vector<double> ran,double r[3]) const
{
  r[0] = x0 + (ix[i] + ran[0])*dx;
  r[1] = y0 + (iy[i] + ran[1])*dy;
  r[2] = z0 + (iz[i] + ran[2])*dz;
}


//------------------------------------------------------------
// return length of zone
//------------------------------------------------------------
double  grid_3D_cart::zone_min_length(const int i) const
{
  return min_ds;
}



//------------------------------------------------------------
// get the velocity vector 
//------------------------------------------------------------
void grid_3D_cart::velocity_vector(const int i, const double x[3], double v[3]) const
{
  // may want to interpolate here
  v[0] = z[i].v[0];
  v[1] = z[i].v[1];
  v[2] = z[i].v[2];
}

//------------------------------------------------------------
// cell-centered coordinates of zone i
//------------------------------------------------------------
void grid_3D_cart::coordinates(const int i,double r[3]) const
{
  r[0] = x0 + (ix[i]+0.5)*dx;
  r[1] = y0 + (iy[i]+0.5)*dy;
  r[2] = z0 + (iz[i]+0.5)*dz;
}


//------------------------------------------------------------
// Write the grid information out to a file
//------------------------------------------------------------
void grid_3D_cart::write_rays(const int iw) const
{
  char zonefile[1000];
  char base[1000];
  int ind;
  double T_rad;
  double r[3];
  ofstream outf;
  int i,j,k;

  if (iw < 10) sprintf(base,"_0000%d",iw);
  else if (iw < 100) sprintf(base,"_000%d",iw);
  else if (iw < 1000) sprintf(base,"_00%d",iw);
  else if (iw < 10000) sprintf(base,"_0%d",iw);
  else sprintf(base,"_%d",iw);

 
  // X-direction
  sprintf(zonefile,"ray%s_x",base);
  outf.open(zonefile);
  outf << setprecision(4);
  outf << scientific;
  outf << "# r[0] r[1] r[2] rho T_rad T_gas Ye" << endl;
  j = ny/2;
  k = nz/2;
  for (i=0;i<nx;i++)
  {
    ind = i*ny*nz + j*nz + k;
    coordinates(ind,r); 
    outf << r[0] << " ";

    T_rad = pow(z[ind].e_rad/pc::a,0.25);
    outf << z[ind].rho   << " ";
    outf << T_rad        << " ";
    outf << z[ind].T_gas << " ";
    outf << z[ind].Ye    << " ";
    outf << endl;
  }
  outf.close();

  // Y-direction
  sprintf(zonefile,"ray%s_y",base);
  outf.open(zonefile);
  outf << setprecision(4);
  outf << scientific;
  outf << "# r[0] r[1] r[2] rho T_rad T_gas Ye" << endl;
  i = nx/2;
  k = nz/2;
  for (j=0; j<ny; j++)
  {
    ind = i*ny*nz + j*nz + k;
    coordinates(ind,r); 
    outf << r[1] << " ";

    T_rad = pow(z[ind].e_rad/pc::a,0.25);
    outf << z[ind].rho   << " ";
    outf << T_rad        << " ";
    outf << z[ind].T_gas << " ";
    outf << z[ind].Ye    << " ";
    outf << endl;
  }
  outf.close();

  // Z-direction
  sprintf(zonefile,"ray%s_z",base);
  outf.open(zonefile);
  outf << setprecision(4);
  outf << scientific;
  outf << "# r[0] r[1] r[2] rho T_rad T_gas Ye" << endl;
  i = nx/2;
  j = ny/2;
  for (k=0; k<nz; k++)
  {
    ind = i*ny*nz + j*nz + k;
    coordinates(ind,r); 
    outf << r[2] << " ";

    T_rad = pow(z[ind].e_rad/pc::a,0.25);
    outf << z[ind].rho   << " ";
    outf << T_rad        << " ";
    outf << z[ind].T_gas << " ";
    outf << z[ind].Ye    << " ";
    outf << endl;
  }
  outf.close();
}


//------------------------------------------------------------
// Reflect off the outer boundary
//------------------------------------------------------------
void grid_3D_cart::reflect_outer(particle *p) const{
  // assumes particle is placed OUTSIDE of the zones

  // invert the radial component of the velocity, put the particle just inside the boundary
  if(p->x[0] < x0){
    assert(p->D[0]<0);
    p->D[0] = -p->D[0];
    p->x[0] = x0 + tiny*dx;
  }
  if(p->x[1] < y0){
    assert(p->D[1]<0);
    p->D[1] = -p->D[1];
    p->x[2] = y0 + tiny*dy;
  }
  if(p->x[2] < z0){
    assert(p->D[0]<0);
    p->D[0] = -p->D[0];
    p->x[1] = z0 + tiny*dz;
  }
  if(p->x[0] > x0+nx*dx){
    assert(p->D[0]>0);
    p->D[0] = -p->D[0];
    p->x[0] = (x0+nx*dx) - tiny*dx;
  }
  if(p->x[1] > y0+ny*dy){
    assert(p->D[1]>0);
    p->D[1] = -p->D[1];
    p->x[1] = (y0+ny*dy) - tiny*dy;
  }
  if(p->x[2] > z0+nz*dz){
    assert(p->D[2]>0);
    p->D[2] = -p->D[2];
    p->x[2] = (z0+nz*dz) - tiny*dz;
  }

  // double check that the particle is in the boundary
  assert(p->x[0]>x0 && p->x[0]<x0+nx*dx);
  assert(p->x[1]>y0 && p->x[1]<y0+ny*dy);
  assert(p->x[2]>z0 && p->x[2]<z0+nz*dz);
}


//------------------------------------------------------------
// Find distance to outer boundary
//------------------------------------------------------------
double grid_3D_cart::dist_to_boundary(const particle *p) const{
  // not yet implemented
  return 0;
}
