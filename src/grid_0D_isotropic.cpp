#include <algorithm>
#include <math.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include "grid_0D_isotropic.h"
#include "global_options.h"

//------------------------------------------------------------
// initialize the zone geometry from model file
//------------------------------------------------------------
void grid_0D_isotropic::read_model_file(Lua* lua)
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
	z.resize(1,zone(3));
	infile >> z[0].rho;
	infile >> z[0].T;
	infile >> z[0].Ye;
	z[0].H_com = 0;
	z[0].e_rad = 0;
	z[0].v[0] = 0;
	z[0].v[1] = 0;
	z[0].v[2] = 0;
	assert(z[0].rho >= 0);
	assert(z[0].T >= 0);
	assert(z[0].Ye >= 0);
	assert(z[0].Ye <= 1.0);
	assert(z[0].v.size() == 3);

	infile.close();
}


//------------------------------------------------------------
// Write a custom model here if you like
//------------------------------------------------------------
void grid_0D_isotropic::custom_model(Lua* lua)
{
	cout << "Error: there is no custom model programmed for grid_0D_isotropic." << endl;
	exit(11);
}


//------------------------------------------------------------
// Return the zone index containing the position x
//------------------------------------------------------------
int grid_0D_isotropic::zone_index(const vector<double>& x) const
{
	return 0;
}


//------------------------------------------------------------
// return volume of zone z_ind
//------------------------------------------------------------
double  grid_0D_isotropic::zone_lab_volume(const int z_ind) const
{
	assert(z_ind == 0);
	return 1.0;
}


//------------------------------------------------------------
// return length of zone
//------------------------------------------------------------
double  grid_0D_isotropic::zone_min_length(const int z_ind) const
{
	assert(z_ind == 0);
	return 1.0;
}


// ------------------------------------------------------------
// find the coordinates of the zone in geometrical coordinates
// ------------------------------------------------------------
void grid_0D_isotropic::zone_coordinates(const int z_ind, vector<double>& r) const{
	assert(z_ind == 0);
	r.resize(dimensionality);
}


//-------------------------------------------
// get directional indices from zone index
//-------------------------------------------
void grid_0D_isotropic::zone_directional_indices(const int z_ind, vector<int>& dir_ind) const
{
	assert(z_ind == 0);
	dir_ind.resize(dimensionality);
}


//------------------------------------------------------------
// sample a random position within the spherical shell
//------------------------------------------------------------
void grid_0D_isotropic::cartesian_sample_in_zone
(const int z_ind, const vector<double>& rand, vector<double>& x) const
{
	assert(z_ind == 0);
	x.resize(3);

	// set the double 3-d coordinates
	x[0] = 0;
	x[1] = 0;
	x[2] = 0;
}


//------------------------------------------------------------
// get the velocity vector 
//------------------------------------------------------------
void grid_0D_isotropic::cartesian_velocity_vector(const vector<double>& x, vector<double>& v) const
{
	assert(x.size()==3);
	v.resize(3);
	v.assign(z[0].v.begin(),z[0].v.end());
	assert(v[0]*v[0] + v[1]*v[1] + v[2]*v[2] <= pc::c*pc::c);
}


//------------------------------------------------------------
// Write the grid information out to a file
//------------------------------------------------------------
void grid_0D_isotropic::write_rays(const int iw) const
{
	// this is a 0D grid, so the function is exactly the same
	// as write_zones
}


//------------------------------------------------------------
// Reflect off the outer boundary
//------------------------------------------------------------
void grid_0D_isotropic::reflect_outer(particle *p) const{
	// does nothing - no boundary
}

//------------------------------------------------------------
// Find distance to outer boundary (less a tiny bit)
// negative distance means inner boundary
//------------------------------------------------------------
double grid_0D_isotropic::lab_dist_to_boundary(const particle *p) const{
	return INFINITY;
}

double grid_0D_isotropic::zone_radius(const int z_ind) const{
	return 0;
}
