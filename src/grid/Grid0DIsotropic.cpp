/*
//  Copyright (c) 2015, California Institute of Technology and the Regents
//  of the University of California, based on research sponsored by the
//  United States Department of Energy. All rights reserved.
//
//  This file is part of Sedonu.
//
//  Sedonu is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  Neither the name of the California Institute of Technology (Caltech)
//  nor the University of California nor the names of its contributors 
//  may be used to endorse or promote products derived from this software
//  without specific prior written permission.
//
//  Sedonu is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with Sedonu.  If not, see <http://www.gnu.org/licenses/>.
//
*/

#include "global_options.h"
#include "Grid0DIsotropic.h"
#include "Lua.h"
#include "Transport.h"

using namespace std;
namespace pc = physical_constants;

Grid0DIsotropic::Grid0DIsotropic(){
	PRINT_ASSERT(NDIMS,==,0);
	grid_type = "Grid0DIsotropic";
}

//------------------------------------------------------------
// initialize the zone geometry from model file
//------------------------------------------------------------
void Grid0DIsotropic::read_model_file(Lua* lua)
{
	// number of zones
	z.resize(1);
	z[0].rho = lua->scalar<double>("Grid0DIsotropic_rho");
	z[0].T   = lua->scalar<double>("Grid0DIsotropic_T")/pc::k_MeV;
	z[0].Ye  = lua->scalar<double>("Grid0DIsotropic_Ye");
	z[0].H_vis = 0;
	PRINT_ASSERT(z[0].rho,>=,0);
	PRINT_ASSERT(z[0].T,>=,0);
	PRINT_ASSERT(z[0].Ye,>=,0);
	PRINT_ASSERT(z[0].Ye,<=,1.0);
}


//------------------------------------------------------------
// Return the zone index containing the position x
//------------------------------------------------------------
int Grid0DIsotropic::zone_index(const double x[3]) const
{
	return 0;
}


//------------------------------------------------------------
// return volume of zone z_ind
//------------------------------------------------------------
double  Grid0DIsotropic::zone_lab_3volume(const int z_ind) const
{
	PRINT_ASSERT(z_ind,==,0);
	return 1.0;
}


//------------------------------------------------------------
// return length of zone
//------------------------------------------------------------
double  Grid0DIsotropic::zone_min_length(const int z_ind) const
{
	PRINT_ASSERT(z_ind,==,0);
	return INFINITY;
}


// ------------------------------------------------------------
// find the coordinates of the zone in geometrical coordinates
// ------------------------------------------------------------
void Grid0DIsotropic::zone_coordinates(const int z_ind, double r[0], const int rsize) const{
	PRINT_ASSERT(z_ind,==,0);
	PRINT_ASSERT(rsize,==,0);
}


//-------------------------------------------
// get directional indices from zone index
//-------------------------------------------
void Grid0DIsotropic::zone_directional_indices(const int z_ind, vector<unsigned>& dir_ind) const
{
	PRINT_ASSERT(z_ind,==,0);
	PRINT_ASSERT(dir_ind.size(),==,0);
}


//------------------------------------------------------------
// sample a random position within the spherical shell
//------------------------------------------------------------
void Grid0DIsotropic::sample_in_zone(const int z_ind, ThreadRNG* rangen, double x[3]) const
{
	PRINT_ASSERT(z_ind,==,0);

	// set the double 3-d coordinates
	x[0] = 0;
	x[1] = 0;
	x[2] = 0;
}


//------------------------------------------------------------
// get the velocity vector 
//------------------------------------------------------------
void Grid0DIsotropic::interpolate_fluid_velocity(const double x[3], double v[3], int z_ind) const
{
	v[0] = 0;
	v[1] = 0;
	v[2] = 0;
}


//------------------------------------------------------------
// Reflect off symmetry plane
//------------------------------------------------------------

void Grid0DIsotropic::symmetry_boundaries(EinsteinHelper *eh, const double tolerance) const{
	// does nothing - no boundary
}

double Grid0DIsotropic::zone_radius(const int z_ind) const{
	return 0;
}

//-----------------------------
// Dimensions of the grid
//-----------------------------
void Grid0DIsotropic::dims(hsize_t dims[0], const int size) const{
	PRINT_ASSERT(size,==,(int)dimensionality());
}

//----------------------------------------------------
// Write the coordinates of the grid points to the hdf5 file
//----------------------------------------------------
void Grid0DIsotropic::write_hdf5_coordinates(H5::H5File file) const
{
	// it's stupid to output this in hdf5...
	cout << "ERROR: write_hdf5_coordinates is not implemented for grid_0D_isotropic." << endl;
	assert(0);
}

double Grid0DIsotropic::lapse(const double xup[4], int z_ind) const {
	return 1.0;
}
void Grid0DIsotropic::shiftup(double betaup[4], const double xup[4], int z_ind) const {
	cout << "g_down not implemented" << endl;
	assert(0);
}
void Grid0DIsotropic::g3_down(const double xup[4], double gproj[4][4], int z_ind) const {
	cout << "g_down not implemented" << endl;
	assert(0);
}
void Grid0DIsotropic::connection_coefficients(const double xup[4], double gamma[4][4][4], int z_ind) const {
	cout << "g_down not implemented" << endl;
	assert(0);
}
double Grid0DIsotropic::zone_lapse(const int z_ind) const{
	return 1.0;
}
void Grid0DIsotropic::axis_vector(vector<Axis>& axes) const{
	axes.resize(0);
}
double Grid0DIsotropic::zone_lorentz_factor(const int z_ind) const{
	return 1.0;
}
