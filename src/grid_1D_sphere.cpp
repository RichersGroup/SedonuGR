#include <algorithm>
#include <math.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include "grid_1D_sphere.h"
#include "global_options.h"

//------------------------------------------------------------
// initialize the zone geometry from model file
//------------------------------------------------------------
void grid_1D_sphere::read_model_file(Lua* lua)
{
	// verbocity
	int my_rank;
	MPI_Comm_rank( MPI_COMM_WORLD, &my_rank );
	const int rank0 = (my_rank == 0);
	if(rank0) cout << "#   Reading 1D model file" << endl;

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
	assert(n_zones > 0);
	z.resize(n_zones,zone(dimensionality));
	r_out.resize(n_zones);

	// read zone properties
	infile >> r_out.min;
	assert(r_out.min >= 0);
	for(int z_ind=0; z_ind<n_zones; z_ind++)
	{
		infile >> r_out[z_ind];
		infile >> z[z_ind].rho;
		infile >> z[z_ind].T;
		infile >> z[z_ind].Ye;
		z[z_ind].H_com = 0;
		z[z_ind].e_rad = 0;
		z[z_ind].v[0] = 0;
		assert(r_out[z_ind] > (z_ind==0 ? r_out.min : r_out[z_ind-1]));
		assert(z[z_ind].rho >= 0);
		assert(z[z_ind].T >= 0);
		assert(z[z_ind].Ye >= 0);
		assert(z[z_ind].Ye <= 1.0);
		assert(z[z_ind].v.size() == dimensionality);
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
// Return the zone index containing the position x
//------------------------------------------------------------
int grid_1D_sphere::zone_index(const vector<double>& x) const
{
	assert(x.size()==3);
	double r = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);

	// check if off the boundaries
	if(r < r_out.min             ) return -1;
	if(r >= r_out[r_out.size()-1] ) return -2;

	// find in zone array using stl algorithm upper_bound and subtracting iterators
	int z_ind = r_out.locate(r);
	assert(z_ind >= 0);
	assert(z_ind < (int)z.size());
	return z_ind;
}


//------------------------------------------------------------
// return volume of zone z_ind
//------------------------------------------------------------
double  grid_1D_sphere::zone_lab_volume(const int z_ind) const
{
	assert(z_ind >= 0);
	assert(z_ind < (int)z.size());
	double r0 = (z_ind==0 ? r_out.min : r_out[z_ind-1]);
	double vol = 4.0*pc::pi/3.0*(r_out[z_ind]*r_out[z_ind]*r_out[z_ind] - r0*r0*r0);
	assert(vol >= 0);
	return vol;
}


//------------------------------------------------------------
// return length of zone
//------------------------------------------------------------
double  grid_1D_sphere::zone_min_length(const int z_ind) const
{
	assert(z_ind >= 0);
	assert(z_ind < (int)z.size());
	double r0 = (z_ind==0 ? r_out.min : r_out[z_ind-1]);
	double min_len = r_out[z_ind] - r0;
	assert(min_len >= 0);
	return min_len;
}


// ------------------------------------------------------------
// find the coordinates of the zone in geometrical coordinates
// ------------------------------------------------------------
void grid_1D_sphere::zone_coordinates(const int z_ind, vector<double>& r) const{
	assert(z_ind >= 0);
	assert(z_ind < (int)z.size());
	r.resize(dimensionality);
	r[0] = 0.5*(r_out[z_ind]+r_out.bottom(z_ind));
	assert(r[0] > 0);
	assert(r[0] < r_out[r_out.size()-1]);
}


//-------------------------------------------
// get directional indices from zone index
//-------------------------------------------
void grid_1D_sphere::zone_directional_indices(const int z_ind, vector<int>& dir_ind) const
{
	assert(z_ind >= 0);
	assert(z_ind < (int)z.size());
	dir_ind.resize(dimensionality);
	dir_ind[0] = z_ind;
}


//------------------------------------------------------------
// sample a random position within the spherical shell
//------------------------------------------------------------
void grid_1D_sphere::cartesian_sample_in_zone
(const int z_ind, const vector<double>& rand, vector<double>& x) const
{
	assert(z_ind >= 0);
	assert(z_ind < (int)z.size());
	assert(rand.size() == 3);
	x.resize(3);

	// inner and outer radii of shell
	double r0 = (z_ind==0 ? r_out.min : r_out[z_ind-1]);
	double r1 = r_out[z_ind];

	// sample radial position in shell using a probability integral transform
	double radius = pow( rand[0]*(r1*r1*r1 - r0*r0*r0) + r0*r0*r0, 1./3.);
	assert(radius >= r0*(1.-tiny));
	assert(radius <= r1*(1.+tiny));
	if(radius<r0) radius = r0;
	if(radius>r1) radius = r1;

	// random spatial angles
	double mu  = 1 - 2.0*rand[1];
	double phi = 2.0*pc::pi*rand[2];
	double sin_theta = sqrt(1 - mu*mu);

	// set the double 3-d coordinates
	x[0] = radius*sin_theta*cos(phi);
	x[1] = radius*sin_theta*sin(phi);
	x[2] = radius*mu;
}


//------------------------------------------------------------
// get the velocity vector 
//------------------------------------------------------------
void grid_1D_sphere::cartesian_velocity_vector(const vector<double>& x, vector<double>& v, int z_ind) const
{
	assert(x.size()==3);
	v.resize(3);
	if(z_ind < 0) z_ind = zone_index(x);
	assert(z_ind >= 0);
	assert(z_ind < (int)z.size());

	// radius in zone
	double r = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);

	// assuming radial velocity (may want to interpolate here)
	// (the other two components are ignored and mean nothing)
	assert(z[z_ind].v.size()==dimensionality);
	v[0] = x[0]/r*z[z_ind].v[0];
	v[1] = x[1]/r*z[z_ind].v[0];
	v[2] = x[2]/r*z[z_ind].v[0];

	// check for pathological case
	if (r == 0)
	{
		v[0] = 0;
		v[1] = 0;
		v[2] = 0;
	}

	assert(v[0]*v[0] + v[1]*v[1] + v[2]*v[2] <= pc::c*pc::c);
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
	double r0 = (r_out.size()>1 ? r_out[r_out.size()-2] : r_out.min);
	double rmax = r_out[r_out.size()-1];
	double dr = rmax - r0;
	double velDotRhat = p->mu();
	double R = p->r();
	assert( fabs(R - r_out[r_out.size()-1]) < tiny*dr);

	// invert the radial component of the velocity
	p->D[0] -= 2.*velDotRhat * p->x[0]/R;
	p->D[1] -= 2.*velDotRhat * p->x[1]/R;
	p->D[2] -= 2.*velDotRhat * p->x[2]/R;
	p->normalize_direction();

	// put the particle just inside the boundary
	double newR = rmax - tiny*dr;
	p->x[0] = p->x[0]/R*newR;
	p->x[1] = p->x[1]/R*newR;
	p->x[2] = p->x[2]/R*newR;

	// must be inside the boundary, or will get flagged as escaped
	assert(zone_index(p->x) >= 0);
}

//------------------------------------------------------------
// Find distance to outer boundary (less a tiny bit)
// negative distance means inner boundary
//------------------------------------------------------------
double grid_1D_sphere::lab_dist_to_boundary(const particle *p) const{
	// Theta = angle between radius vector and direction (Pi if outgoing)
	// Phi   = Pi - Theta (angle on the triangle) (0 if outgoing)
	double Rout  = r_out[r_out.size()-1];
	double Rin   = r_out.min;
	double r  = p->r();
	double mu = p->mu();
	double d_outer_boundary = numeric_limits<double>::infinity();
	double d_inner_boundary = numeric_limits<double>::infinity();
	assert(r<Rout);
	assert(zone_index(p->x) >= -1);

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
	if(d_inner_boundary<=0 && fabs(d_inner_boundary/Rin)<tiny*(r_out[0]-Rin)) d_inner_boundary = tiny*(r_out[0]-Rin);
	assert(d_inner_boundary > 0);

	// distance to outer boundary
	d_outer_boundary = -r*mu + sqrt(r*r*(mu*mu-1.0) + Rout*Rout);
	if(d_outer_boundary<=0 && fabs(d_outer_boundary/Rin)<tiny*(Rout-r_out[r_out.size()-1])) d_outer_boundary = tiny*(Rout-r_out[r_out.size()-1]);
	assert(d_outer_boundary > 0);
	assert(d_outer_boundary <= 2.*Rout);

	// make sure the particle ends up in a reasonable place
	return min(d_inner_boundary, d_outer_boundary);
}

double grid_1D_sphere::zone_radius(const int z_ind) const{
	assert(z_ind >= 0);
	assert(z_ind < (int)z.size());
	return r_out[z_ind];
}
