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

#include <mpi.h>
#include <fstream>
#include <sstream>
#include "Grid3DCart.h"
#include "global_options.h"
#include "Transport.h"
#include "H5Cpp.h"

using namespace std;
namespace pc = physical_constants;

//------------
// constructor
//------------
Grid3DCart::Grid3DCart(){
	for(int i=0; i<3; i++){
		nx[i]      = -1;
		dx[i]      = NaN;
		x0[i]      = NaN;
		x1[i]      = NaN;
		xmax[i]    = NaN;
		reflect[i] = 0;
	}
	rotate_hemisphere[0] = 0;
	rotate_hemisphere[1] = 0;
	rotate_quadrant = 0;
}

//------------------------------------------------------------
// Read in a cartesian model file
//------------------------------------------------------------
void Grid3DCart::read_model_file(Lua* lua)
{
	std::string model_type = lua->scalar<std::string>("model_type");
	if(model_type == "SpEC") read_SpEC_file(lua);
	else if(model_type == "THC") read_THC_file(lua);
	else{
		cout << "ERROR: model type unknown." << endl;
		exit(8);
	}
}



void Grid3DCart::read_THC_file(Lua* lua)
{
	// verbocity
	int my_rank;
	MPI_Comm_rank( MPI_COMM_WORLD, &my_rank );
	const int rank0 = (my_rank == 0);

	// conversion factors
	double convert_length = 147690.2071535873; // code --> cm
	double convert_density = 6.173937319029555e+17; // code --> g/ccm
	double convert_velocity = pc::c; // code --> cm/s
	double convert_temperature = 1./pc::k_MeV; // code --> K

	// generalHDF5 variables
	H5::DataSet dataset;
	H5::DataSpace space;
	H5::CompType comptype;
	H5::Attribute attr;
	H5::Group group;

	// open the model files
	if(rank0) cout << "# Reading the model file..." << endl;
	string model_filename   = lua->scalar<string>("model_file"  );
	H5::H5File file(model_filename, H5F_ACC_RDONLY);

	// get the refinement level
	int reflevel = lua->scalar<int>("Grid3DCart_THC_reflevel");
	stringstream groupname;
	groupname << "/reflevel=" << reflevel << "/";
	group = file.openGroup(groupname.str());

	// open Ye
	dataset = file.openDataSet(groupname.str()+"Ye");
	space = dataset.getSpace();
	comptype = dataset.getCompType();

	// get the dataset dimensions
	const int dataset_rank = 3;
	PRINT_ASSERT(space.getSimpleExtentNdims(),==,dataset_rank);                 // 3D array
	hsize_t hdf5_dims[dataset_rank];
	space.getSimpleExtentDims(hdf5_dims);
	for(int i=0; i<3; i++) nx[i] = hdf5_dims[i];
	int dataset_nzones = nx[0] * nx[1] * nx[2];

	//======================//
	// read the coordinates //
	//======================//
	// read in the deltas
	attr = group.openAttribute("delta");
	attr.read(H5::PredType::IEEE_F64LE,&dx[0]);
	for(int i=0; i<3; i++){
		dx[i] *= convert_length;
		PRINT_ASSERT(dx[i],>,0);
	}

	// read in the extent.
	double extent[6];
	attr = group.openAttribute("extent");
	attr.read(H5::PredType::IEEE_F64LE,extent);
	for(int i=0; i<6; i++) extent[i] *= convert_length;

	// set the grid structure variables. half-cell offset since data is vertex-centered.
	for(int i=0; i<3; i++){
		x0[i] = extent[2*i] - dx[i]/2.0;
		x1[i] = x0[i] + dx[i];
		xmax[i] = extent[2*i+1] + dx[i]/2.0;
		PRINT_ASSERT(xmax[i],>,x0[i]);

		// check that the number of data points is consistent with the range and delta
		double ntmp = (xmax[i] - x0[i])/dx[i];
		PRINT_ASSERT(abs(nx[i]-ntmp),<,TINY);
	}

	// modify grid in event of symmetry
	reflect[0] = lua->scalar<int>("Grid3DCart_reflect_x");
	reflect[1] = lua->scalar<int>("Grid3DCart_reflect_y");
	reflect[2] = lua->scalar<int>("Grid3DCart_reflect_z");
	rotate_hemisphere[0] = lua->scalar<int>("Grid3DCart_rotate_hemisphere_x"); // rotate around z
	rotate_hemisphere[1] = lua->scalar<int>("Grid3DCart_rotate_hemisphere_y"); // rotate around z
	rotate_quadrant = lua->scalar<int>("Grid3DCart_rotate_quadrant"); // rotate around z, x and y both positive

	// check parameters
	if(reflect[0] || reflect[1]){
		assert(!rotate_quadrant);
		assert(!rotate_hemisphere[0]);
		assert(!rotate_hemisphere[1]);
	}
	if(rotate_quadrant){
		assert(!rotate_hemisphere[0]);
		assert(!rotate_hemisphere[1]);
	}
	assert(!( rotate_hemisphere[0] && rotate_hemisphere[1] ));

	// decide which axes to truncate
	bool truncate[3] = {false,false,false};
	int offset[3] = {0,0,0};
	truncate[0] = reflect[0] || rotate_hemisphere[0] || rotate_quadrant;
	truncate[1] = reflect[1] || rotate_hemisphere[1] || rotate_quadrant;
	truncate[2] = reflect[2];

	// truncate the grid to create appropriate boundaries
	for(int i=0; i<3; i++) if(truncate[i]){
		PRINT_ASSERT(xmax[i],>,0);
		PRINT_ASSERT(x0[i],<=,0);
		x0[i] = 0.0;
		x1[i] = fmod(xmax[i],dx[i]);
		if(abs(x1[i]-x0[i])/dx[i] < TINY) x1[i] += dx[i]; // don't let the leftmost grid cell be stupidly TINY.
		PRINT_ASSERT(fmod(xmax[i]-x1[i],dx[i]) < TINY,||,fmod(xmax[i]-x1[i],dx[i])-dx[i] < TINY); // make sure cells line up
		nx[i] = (int)((xmax[i]-x1[i])/dx[i] + 0.5) +1; // 0.5 to deal with finite precision.
		offset[i] = hdf5_dims[i] - nx[i];
	}
	// truncate outer part of zones if quadrant rotational symmetry
	if(rotate_quadrant) for(int i=0; i<2; i++){
		int iother = (i+1)%2;
		if(nx[i] > nx[iother]){
			PRINT_ASSERT(x0[i]-x0[iother],==,0);
			PRINT_ASSERT(x1[i]-x1[iother],<,TINY);
			PRINT_ASSERT(dx[i]-dx[iother],<,TINY);
			nx[i] = nx[iother];
			xmax[i] = xmax[iother];
		}
	}	
	
	// set up the zone structure
	int nzones = nx[0] * nx[1] * nx[2];
	z.resize(nzones);

	// print out grid structure
	if(rank0){
		cout << "#   Using refinement level " << reflevel << endl;
		cout << "#   Minima : {" <<   x0[0] << "," <<   x0[1] << "," <<   x0[2] << "} cm" << endl;
		cout << "#   Next   : {" <<   x1[0] << "," <<   x1[1] << "," <<   x1[2] << "} cm" << endl;
		cout << "#   Maxima : {" << xmax[0] << "," << xmax[1] << "," << xmax[2] << "} cm" << endl;
		cout << "#   Deltas : {" <<   dx[0] << "," <<   dx[1] << "," <<   dx[2] << "} cm" << endl;
		cout << "#   Number : {" <<   nx[0] << "," <<   nx[1] << "," <<   nx[2] << "}"    << endl;
	}


	//=========================//
	// read in the actual data //
	//=========================//
	vector<double>   Ye(dataset_nzones,0.0);
	vector<double>  rho(dataset_nzones,0.0);
	vector<double> temp(dataset_nzones,0.0);
	vector<double> velx(dataset_nzones,0.0);
	vector<double> vely(dataset_nzones,0.0);
	vector<double> velz(dataset_nzones,0.0);
	vector<double>  vol(dataset_nzones,0.0);
	dataset = file.openDataSet(groupname.str()+"Ye");
	dataset.read(&Ye[0],H5::PredType::IEEE_F64LE);
	dataset = file.openDataSet(groupname.str()+"rho");
	dataset.read(&(rho[0]),H5::PredType::IEEE_F64LE);
	dataset = file.openDataSet(groupname.str()+"temp");
	dataset.read(&(temp[0]),H5::PredType::IEEE_F64LE);
	dataset = file.openDataSet(groupname.str()+"velx");
	dataset.read(&(velx[0]),H5::PredType::IEEE_F64LE);
	dataset = file.openDataSet(groupname.str()+"vely");
	dataset.read(&(vely[0]),H5::PredType::IEEE_F64LE);
	dataset = file.openDataSet(groupname.str()+"velz");
	dataset.read(&(velz[0]),H5::PredType::IEEE_F64LE);
	dataset = file.openDataSet(groupname.str()+"vol");
	dataset.read(&(vol[0]),H5::PredType::IEEE_F64LE);
	#pragma omp parallel for
	for(int i=0; i<dataset_nzones; i++){
		rho[i] *= convert_density;
		temp[i] *= convert_temperature;
		velx[i] *= convert_velocity;
		vely[i] *= convert_velocity;
		velz[i] *= convert_velocity;
	}


	//===============//
	// fill the grid //
	//===============//
    #pragma omp parallel for
	for(int z_ind=0; z_ind<z.size(); z_ind++){

		// directional indices in Sedonu grid
		int dir_ind[3];
		zone_directional_indices(z_ind,dir_ind,3);

		// directional indices in hdf5 data
		unsigned hdf5_dir_ind[3];
		for(int d=0; d<3; d++){
			hdf5_dir_ind[d] = offset[d] + dir_ind[d];
			PRINT_ASSERT(hdf5_dir_ind[d],>=,0);
			PRINT_ASSERT(offset[d],>=,0);
		}

		// global hdf5 index
		const int dataset_ind = hdf5_dir_ind[2] + hdf5_dims[2]*hdf5_dir_ind[1] + hdf5_dims[1]*hdf5_dims[2]*hdf5_dir_ind[0];
		PRINT_ASSERT((int)dataset_ind,<,(int)(hdf5_dims[0] * hdf5_dims[1] * hdf5_dims[2]));
		PRINT_ASSERT((int)dataset_ind,>=,0);

		// fill the zone
		z[z_ind].rho  =  rho[dataset_ind];
		z[z_ind].T    = temp[dataset_ind];
		z[z_ind].Ye   =   Ye[dataset_ind];
		z[z_ind].u[0] = velx[dataset_ind];
		z[z_ind].u[1] = vely[dataset_ind];
		z[z_ind].u[2] = velz[dataset_ind];

		PRINT_ASSERT(z[z_ind].rho,>=,0.0);
		PRINT_ASSERT(z[z_ind].T,>=,0.0);
		PRINT_ASSERT(z[z_ind].Ye,>=,0.0);
		PRINT_ASSERT(z[z_ind].Ye,<=,1.0);
		PRINT_ASSERT(zone_speed2(z_ind),<,pc::c*pc::c);
		PRINT_ASSERT(zone_speed2(z_ind),>=,0.0);
	}

	file.close();
}

void Grid3DCart::read_SpEC_file(Lua* lua)
{
	// get mpi rank
	int my_rank;
	MPI_Comm_rank( MPI_COMM_WORLD, &my_rank  );

	// remember which axes were truncated
	vector<bool> truncate = vector<bool>(3,false);
	truncate[0] = reflect[0] || rotate_hemisphere[0] || rotate_quadrant;
	truncate[1] = reflect[1] || rotate_hemisphere[1] || rotate_quadrant;
	truncate[2] = reflect[2];

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
	for(int i=0; i<3; i++) infile >> nx[i];
	// assume reflecting symmetry?
	for(int i=0; i<3; i++){
		infile >> reflect[i];
		if(reflect[i]) nx[i]*=2;
		else if(i<2 && rotate_hemisphere[i]) nx[i]*=2;
		else if(i<2 && rotate_quadrant) nx[i]*=4;
	}
	int n_zones = nx[0]*nx[1]*nx[2];

	//set the zone sizes and volumes
	for(int i=0; i<3; i++){
		infile >> x0[i];
		infile >> xmax[i];
		dx[i] = (xmax[i]-x0[i])/(double)nx[i];
		if(reflect[i]) dx[i]*=2;
		else if(i<2 && rotate_hemisphere[i]) dx[i]*=2;
		else if(i<2 && rotate_quadrant) dx[i]*=4;
		x1[i] = x0[i]+dx[i];
	}

	// First loop - set indices and read zone values in file
	// loop order is the file order
	// index order matches that in get_zone
	z.resize(n_zones);
	int ind = 0;
	bool rx,ry,rz;
	for (int k=0;k<nx[2];k++)
		for (int j=0;j<nx[1];j++)
			for (int i=0;i<nx[0];i++)
			{
				// set current index
				ind = zone_index(i,j,k);

				// create reverse map to x,y,z indices
				rx=false; ry=false; rz=false;
				if(truncate[0] && i<nx[0]/2) rx = true;
				if(truncate[1] && j<nx[1]/2) ry = true;
				if(truncate[2] && k<nx[2]/2) rz = true;

				// read in values if not in reflected/rotated zone
				if(!rx && !ry && !rz)
				{
					infile >> z[ind].rho;
					infile >> z[ind].T;
					infile >> z[ind].Ye;
					infile >> z[ind].u[0];
					infile >> z[ind].u[1];
					infile >> z[ind].u[2];
				}
				else{ //poison values
					z[ind].rho   = NaN;
					z[ind].T     = NaN;
					z[ind].Ye    = NaN;
					z[ind].u[0]  = NaN;
					z[ind].u[1]  = NaN;
					z[ind].u[2]  = NaN;
				}
			}

	// Second loop - apply symmetries
	ind=0;
	int origin_ind, origin_i, origin_j, origin_k;
	for (int k=0;k<nx[2];k++)
		for (int j=0;j<nx[1];j++)
			for (int i=0;i<nx[0];i++)
			{
				// set current index
				ind = zone_index(i,j,k);

				// are we in a reflection/rotation zone?
				rx=false; ry=false; rz=false;
				if(truncate[0] && i<nx[0]/2) rx = 1;
				if(truncate[1] && j<nx[1]/2) ry = 1;
				if(truncate[2] && k<nx[2]/2) rz = 1;

				// copy appropriate values
				if(rx) origin_i = (nx[0]-1)-i; else origin_i = i;
				if(ry) origin_j = (nx[1]-1)-j; else origin_j = j;
				if(rz) origin_k = (nx[2]-1)-k; else origin_k = k;
				if(rx || ry || rz)
				{
					origin_ind = zone_index(origin_i,origin_j,origin_k);
					z[ind].rho   = z[origin_ind].rho;
					z[ind].T = z[origin_ind].T;
					z[ind].Ye    = z[origin_ind].Ye;
					z[ind].u[0]  = z[origin_ind].u[0];
					z[ind].u[1]  = z[origin_ind].u[1];
					z[ind].u[2]  = z[origin_ind].u[2];
				}
			}

	// adjust x0,y0,z0 to indicate new, reflected lower boundary
	if(my_rank==0) cout << "# (zmin,zmax) before adjusting: (" << x0[2] << "," << xmax[2] << ")" << endl;
	for(int i=0; i<3; i++) if(truncate[i]) x0[i] = x0[i] - (xmax[i]-x0[i]);

	// debugging some output
	if(my_rank==0){
		cout << "#   nx=" << nx[0] << endl << "# ny=" << nx[1] << endl << "# nz=" << nx[2] << endl;
		cout << "#   number of zones:" << z.size() << endl;
		cout << "#   minima:{" << x0[0] << ", " << x0[1] << ", " << x0[2] << "}" << endl;
		cout << "#   maxima:{" << x0[0]+(nx[0]*dx[0]) << ", " << x0[1]+(nx[1]*dx[1]) << ", " << x0[2]+(nx[2]*dx[2]) << "}" << endl;
		cout << "#   deltas:{" << dx[0] << ", " << dx[1] << ", " << dx[2] << "}" << endl;
	}
}


//------------------------------------------------------------
// Overly simple search to find zone
//------------------------------------------------------------
int Grid3DCart::zone_index(const double x[3], const int xsize) const
{
	PRINT_ASSERT(xsize,==,3);

	// check for off grid
	for(int i=0; i<3; i++) if (x[i]<x0[i] || x[i]>xmax[i]) return -2;

	// get directional indices
	int dir_ind[3] = {-1,-1,-1};
	for(int i=0; i<3; i++){
		if(x[i]<=x1[i]) dir_ind[i] = 0;
		else if(x[i] > xmax[i]-dx[i]) dir_ind[i] = nx[i]-1; //need this to prevent issues when particle is ON boundary
		else dir_ind[i] = (x[i]-x1[i]) / dx[i] + 1.0;
		PRINT_ASSERT(dir_ind[i],>=,0);
		PRINT_ASSERT(dir_ind[i],<,nx[i]);
	}

	int z_ind = zone_index(dir_ind[0],dir_ind[1],dir_ind[2]);
	PRINT_ASSERT(z_ind,>=,0);
	PRINT_ASSERT(z_ind,<,(int)z.size());
	return z_ind;
}


//-------------------------------------------------
// get the zone index from the directional indices
//-------------------------------------------------
int Grid3DCart::zone_index(const int i, const int j, const int k) const{
	PRINT_ASSERT(i,>=,0);
	PRINT_ASSERT(j,>=,0);
	PRINT_ASSERT(k,>=,0);
	PRINT_ASSERT(i,<,nx[0]);
	PRINT_ASSERT(j,<,nx[1]);
	PRINT_ASSERT(k,<,nx[2]);
	const int z_ind = i*nx[1]*nx[2] + j*nx[2] + k;
	PRINT_ASSERT(z_ind,<,(int)z.size());
	return z_ind;
}

//-------------------------------------------
// get directional indices from zone index
//-------------------------------------------
void Grid3DCart::zone_directional_indices(const int z_ind, int dir_ind[3], const int size) const
{
	PRINT_ASSERT(z_ind,>=,0);
	PRINT_ASSERT(z_ind,<,(int)z.size());
	PRINT_ASSERT(size,==,dimensionality());

	dir_ind[0] =  z_ind / (nx[1]*nx[2]);
	dir_ind[1] = (z_ind % (nx[1]*nx[2])) / nx[2];
	dir_ind[2] =  z_ind % nx[2];

	for(int i=0; i<3; i++){
		PRINT_ASSERT(dir_ind[i],>=,0);
		PRINT_ASSERT(dir_ind[i],<,nx[i]);
	}
}


//------------------------------------------------------------
// return volume of zone (precomputed)
//------------------------------------------------------------
double Grid3DCart::zone_lab_volume(const int z_ind) const
{
	PRINT_ASSERT(z_ind,>=,0);
	PRINT_ASSERT(z_ind,<,(int)z.size());
	int dir_ind[3];
	zone_directional_indices(z_ind,dir_ind,3);
	double delta[3];
	get_deltas(z_ind,delta,3);
	return delta[0] * delta[1] * delta[2];
}


//------------------------------------------------------------
// sample a random position within the cubical cell
//------------------------------------------------------------
void Grid3DCart::cartesian_sample_in_zone
(const int z_ind, const double rand[3], const int randsize, double x[3], const int xsize) const
{
	PRINT_ASSERT(z_ind,>=,0);
	PRINT_ASSERT(z_ind,<,(int)z.size());
	PRINT_ASSERT(xsize,==,3);

	// zone directional indices
	int dir_ind[3];
	zone_directional_indices(z_ind,dir_ind,3);

	// zone deltas in each of three directions
	double delta[3];
	get_deltas(z_ind,delta,3);

	// set the random location
	for(int i=0; i<3; i++){
		x[i] = zone_left_boundary(i,dir_ind[i]) + delta[i]*rand[i];

		// make sure the particle is in bounds
		PRINT_ASSERT(x[i],>,x0[i]   - TINY*dx[i]);
		PRINT_ASSERT(x[0],<,xmax[i] + TINY*dx[i]);

		// return particles just outside cell boundaries to within cell boundaries
		x[i] = min(x[i], zone_right_boundary(i,dir_ind[i]));
		x[i] = max(x[i],  zone_left_boundary(i,dir_ind[i]));
	}
}


//------------------------------------------------------------
// return length of zone
//------------------------------------------------------------
double  Grid3DCart::zone_min_length(const int z_ind) const
{
	PRINT_ASSERT(z_ind,>=,0);
	PRINT_ASSERT(z_ind,<,(int)z.size());

	double delta[3];
	get_deltas(z_ind,delta,3);

	double min_ds = min(delta[0], min(delta[1],delta[2]) );
	return min_ds;
}



//------------------------------------------------------------
// get the velocity vector 
//------------------------------------------------------------
void Grid3DCart::cartesian_velocity_vector(const double x[3], const int xsize, double v[3], const int vsize, int z_ind) const
{
	PRINT_ASSERT(xsize,==,3);
	PRINT_ASSERT(vsize,==,3);
	if(z_ind<0) z_ind = zone_index(x,xsize);
	PRINT_ASSERT(z_ind,>=,0);
	PRINT_ASSERT(z_ind,<,(int)z.size());

	// may want to interpolate here?
	for(int i=0; i<vsize; i++) v[i] = z[z_ind].u[i];

	PRINT_ASSERT(v[0]*v[0] + v[1]*v[1] + v[2]*v[2],<=,pc::c*pc::c);
}

//------------------------------------------------------------
// cell-centered coordinates of zone i
//------------------------------------------------------------
void Grid3DCart::zone_coordinates(const int z_ind, double r[3], const int rsize) const
{
	PRINT_ASSERT(z_ind,>=,0);
	PRINT_ASSERT(z_ind,<,(int)z.size());
	PRINT_ASSERT(rsize,==,dimensionality());
	int dir_ind[3];
	zone_directional_indices(z_ind,dir_ind,3);

	for(int i=0; i<3; i++)
		r[i] = ( dir_ind[i]==0 ? 0.5*(x0[i]+x1[i]) : x1[i] + ((double)(dir_ind[i]-1) + 0.5) * dx[i] );
}


//------------------------------------------------------------
// Write the grid information out to a file
//------------------------------------------------------------
void Grid3DCart::write_rays(const int iw) const
{

	int i,j,k;
	int z_ind;
	double r[3];
	string filename = "";
	ofstream outf;

	// get the origin coordinates
	double origin[3] = {0,0,0};
	z_ind = zone_index(origin,3);
	int iorigin[3] = {0,0,0};
	zone_directional_indices(z_ind,iorigin,3);


	// XY-slice
	filename = Transport::filename("slice_xy",iw,".dat");
	outf.open(filename.c_str());
	write_header(outf);
	k = iorigin[2]; //nx[2]/2;
	for(i=0; i<nx[0]; i++) for(j=0; j<nx[1]; j++){
		if(j==0) outf << endl;
		z_ind = zone_index(i,j,k);
		zone_coordinates(z_ind,r,3);
		write_line(outf,z_ind);
	}
	outf.close();

	// XZ-slice
	filename = Transport::filename("slice_xz",iw,".dat");
	outf.open(filename.c_str());
	write_header(outf);
	j = iorigin[1]; //nx[1]/2;
	for(i=0; i<nx[0]; i++) for(k=0; k<nx[2]; k++){
		if(k==0) outf << endl;
		z_ind = zone_index(i,j,k);
		zone_coordinates(z_ind,r,3);
		write_line(outf,z_ind);
	}
	outf.close();

	// YZ-slice
	filename = Transport::filename("slice_yz",iw,".dat");
	outf.open(filename.c_str());
	write_header(outf);
	i = iorigin[0]; //nx[0]/2;
	for(j=0; j<nx[1]; j++) for(k=0; k<nx[2]; k++){
		if(k==0) outf << endl;
		z_ind = zone_index(i,j,k);
		zone_coordinates(z_ind,r,3);
		write_line(outf,z_ind);
	}
	outf.close();

	// X-direction
	filename = Transport::filename("ray_x",iw,".dat");
	outf.open(filename.c_str());
	write_header(outf);
	j = nx[1]/2;
	k = nx[2]/2;
	for (i=0;i<nx[0];i++){
		z_ind = zone_index(i,j,k);
		zone_coordinates(z_ind,r,3);
		write_line(outf,z_ind);
	}
	outf.close();

	// Y-direction
	filename = Transport::filename("ray_y",iw,".dat");
	outf.open(filename.c_str());
	write_header(outf);
	i = nx[0]/2;
	k = nx[2]/2;
	for (j=0; j<nx[1]; j++){
		z_ind = zone_index(i,j,k);
		zone_coordinates(z_ind,r,3);
		write_line(outf,z_ind);
	}
	outf.close();

	// Z-direction
	filename = Transport::filename("ray_z",iw,".dat");
	outf.open(filename.c_str());
	write_header(outf);
	i = nx[0]/2;
	j = nx[1]/2;
	for (k=0; k<nx[2]; k++)
	{
		z_ind = zone_index(i,j,k);
		zone_coordinates(z_ind,r,3);
		write_line(outf,z_ind);
	}
	outf.close();
}


//------------------------------------------------------------
// Reflect off the outer boundary
//------------------------------------------------------------
void Grid3DCart::reflect_outer(LorentzHelper *lh) const{
	const Particle *p = lh->particle_readonly(lab);

	// assumes particle is placed OUTSIDE of the zones
	int z_ind = zone_index(p->xup,3);
	double delta[3];
	get_deltas(z_ind,delta,3);

	double Dlab[3];
	lh->p_D(lab,Dlab,3);

	// invert the radial component of the velocity, put the particle just inside the boundary
	double x3[3], D[3];
	for(int i=0; i<3; i++){
		// initialize values
		x3[i] = lh->p_xup()[i];
		D[i] = Dlab[i];

		// inner boundary
		if(p->xup[i] < x0[i]){
			PRINT_ASSERT(Dlab[i],<,0);
			D[i] = -Dlab[i];
			x3[i] = x0[i] + TINY*delta[i];
		}

		// outer boundary
		if(p->xup[i] > xmax[i]){
			PRINT_ASSERT(Dlab[i],>,0);
			D[i] = -Dlab[i];
			x3[i] = xmax[i] - TINY*delta[i];
		}

		// double check that the particle is in the boundary
		PRINT_ASSERT(x3[i],>,x0[i]);
		PRINT_ASSERT(x3[i],<,xmax[i]);

		// assign the arrays
		double xup[4] = {x3[0], x3[1], x3[2], p->xup[3]};
		lh->set_p_xup(xup,4);
		lh->set_p_D<lab>(D,3);
	}
}


//------------------------------------------------------------
// Find distance to outer boundary
//------------------------------------------------------------
double Grid3DCart::lab_dist_to_boundary(const LorentzHelper *lh) const{
	const Particle *p = lh->particle_readonly(lab);
	double Dlab[3];
	lh->p_D(lab,Dlab,3);

	bool inside = true;
	for(int i=0; i<3; i++) inside = inside && (p->xup[i] >= x0[i]) && (p->xup[i] <= xmax[i]);

	// case: particle is inside the boundaries
	if(inside){
		double dist[3] = {-1,-1,-1};
		for(int i=0; i<3; i++){
			dist[i] = (Dlab[i]>0 ? xmax[i]-p->xup[i] : p->xup[i]-x0[i]);
			dist[i] /= fabs(Dlab[i]);
			PRINT_ASSERT(dist[i],>=,0);
		}

		double final_dist = min(min(dist[0],dist[1]),dist[2]);
		PRINT_ASSERT(final_dist,<=,sqrt((xmax[0]-x0[0])*(xmax[0]-x0[0]) + (xmax[1]-x0[1])*(xmax[1]-x0[1]) + (xmax[2]-x0[2])*(xmax[2]-x0[2])) );
		return final_dist;
	}

	// case: particle is outside the boundaries.
	else return numeric_limits<double>::infinity();
}


double Grid3DCart::zone_radius(const int z_ind) const{
	PRINT_ASSERT(z_ind,>=,0);
	PRINT_ASSERT(z_ind,<,(int)z.size());
	double r[3];
	zone_coordinates(z_ind,r,3);
	return sqrt(dot_Minkowski<3>(r,r,3));
}

//-----------------------------
// Dimensions of the grid
//-----------------------------
void Grid3DCart::dims(hsize_t dims[3], const int size) const{
	PRINT_ASSERT(size,==,dimensionality());
	for(int i=0; i<3; i++) dims[i] = nx[i];
}

//----------------------------------------------------
// Write the coordinates of the grid points to the hdf5 file
//----------------------------------------------------
void Grid3DCart::write_hdf5_coordinates(H5::H5File file) const
{
	// useful quantities
	H5::DataSet dataset;
	H5::DataSpace dataspace;
	vector<float> tmp;

	// get dimensions
	hsize_t coord_dims[3];
	dims(coord_dims,3);
	for(unsigned i=0; i<3; i++) coord_dims[i]++; //make room for min value

	// write x coordinates
	for(int dir=0; dir<3; dir++){
		dataspace = H5::DataSpace(1,&coord_dims[dir]);
		stringstream dirstream;
		dirstream << dir;
		dataset = file.createDataSet("grid_"+dirstream.str()+"(cm)",H5::PredType::IEEE_F32LE,dataspace);
		tmp.resize(coord_dims[dir]);
		tmp[0] = x0[dir];
		tmp[1] = x1[dir];
		if(nx[dir]>1) for(int i=2; i<nx[dir]+1; i++) tmp[i] = x1[dir] + (i-1)*dx[dir];
		dataset.write(&tmp[0],H5::PredType::IEEE_F32LE);
		dataset.close();
	}
}


void Grid3DCart::get_deltas(const int z_ind, double delta[3], const int size) const
{
	PRINT_ASSERT(z_ind,<,(int)z.size());
	PRINT_ASSERT(size,==,3);

	// get directional indices
	int dir_ind[3];
	zone_directional_indices(z_ind,dir_ind,3);

	for(int i=0; i<3; i++){
		delta[i] = (dir_ind[i]==0 ? x1[i]-x0[i] : dx[i]);
		PRINT_ASSERT(delta[i],>,0);
	}
}


double Grid3DCart::zone_left_boundary(const unsigned dir, const unsigned dir_ind) const{
	PRINT_ASSERT(dir,>=,0);
	PRINT_ASSERT(dir,<,3);
	PRINT_ASSERT(dir_ind,>=,0);

	double boundary = ( dir_ind==0 ? x0[dir] : x1[dir] + (double)(dir_ind-1) * dx[dir] );
	PRINT_ASSERT(boundary,<=,xmax[dir]);
	PRINT_ASSERT(boundary,>=,x0[dir]);
	return boundary;
}
double Grid3DCart::zone_right_boundary(const unsigned dir, const unsigned dir_ind) const{
	PRINT_ASSERT(dir,>=,0);
	PRINT_ASSERT(dir,<,3);
	PRINT_ASSERT(dir_ind,>=,0);

	double boundary = ( dir_ind==0 ? x1[dir] : x1[dir] + (double)dir_ind * dx[dir] );
	PRINT_ASSERT(boundary,<=,xmax[dir]*(1.0+TINY));
	PRINT_ASSERT(boundary,>=,x0[dir]);
	return boundary;
}


//------------------------------------------------------------
// Reflect off revlecting boundary condition
//------------------------------------------------------------
void Grid3DCart::symmetry_boundaries(LorentzHelper *lh) const{

	double D[3], x3[3];

	// initialize the arrays
	for(int i=0; i<3; i++) x3[i] = lh->p_xup()[i];
	lh->p_D(lab,D,3);

	// invert the radial component of the velocity, put the particle just inside the boundary
	for(int i=0; i<3; i++){
		if(reflect[i] && x3[i] < x0[i]){
			PRINT_ASSERT(x0[i]-x3[i],<,TINY*dx[i]);
			PRINT_ASSERT(x0[i],==,0);
			PRINT_ASSERT(D[i],<,0);
			D[i] = -D[i];
			x3[i] = x0[i] + TINY*(x1[i]-x0[i]);

			// double check that the particle is in the boundary
			PRINT_ASSERT(x3[i],>=,x0[i]);
			PRINT_ASSERT(x3[i],<=,xmax[i]);
		}
	}

	// rotating boundary conditions
	for(int i=0; i<2; i++){
		if(x3[i] < x0[i] && (rotate_hemisphere[i] || rotate_quadrant)){
			PRINT_ASSERT(x0[i]-x3[i],<,TINY*(x1[i]-x0[i]));
			PRINT_ASSERT(x0[i],==,0);
			PRINT_ASSERT(D[i],<,0);
			int other = i==0 ? 1 : 0;
			
			if(rotate_hemisphere[i]){
				for(int j=0; j<2; j++) D[j] = -D[j];
				x3[i    ] = x0[i] + TINY*(x1[i]-x0[i]);
				x3[other] = -x3[other];
			}
			else if(rotate_quadrant){
				double tmp = D[i];
				D[i] = D[other];
				D[other] = -tmp;
				x3[i] = x3[other];
				x3[other] = x0[other] + TINY*(x1[other]-x0[other]);
			}

			// double check that the particle is in the boundary
			for(int j=0; j<3; j++){
				PRINT_ASSERT(x3[j],>=,x0[j]);
				PRINT_ASSERT(x3[j],<=,xmax[j]);
			}
		}

		// assign the arrays
		double xup[4] = {x3[0], x3[1], x3[2], lh->p_xup()[3]};
		lh->set_p_xup(xup,4);
		lh->set_p_D<lab>(D,3);
	}
}
