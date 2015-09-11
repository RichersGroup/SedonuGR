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

#include "grid_3D_cart.h"
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include "global_options.h"
#include "mpi.h"
#include "Lua.h"

//------------
// constructor
//------------
grid_3D_cart::grid_3D_cart(){
	nx=-MAXLIM;ny=-MAXLIM;nz=-MAXLIM;
	dx=NaN;dy=NaN;dz=NaN;
	x0=NaN;y0=NaN;z0=NaN;
	x1=NaN;y1=NaN;z1=NaN;
	xmax=NaN;ymax=NaN;zmax=NaN;
	reflect_x=0;reflect_y=0;reflect_z=0;
}

//------------------------------------------------------------
// Read in a cartesian model file
//------------------------------------------------------------
void grid_3D_cart::read_model_file(Lua* lua)
{
	std::string model_type = lua->scalar<std::string>("model_type");
	if(model_type == "SpEC") read_SpEC_file(lua);
	else if(model_type == "David") read_David_file(lua);
	else{
		cout << "ERROR: model type unknown." << endl;
		exit(8);
	}
}



void grid_3D_cart::read_David_file(Lua* lua)
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
	int reflevel = lua->scalar<int>("reflevel");
	stringstream groupname;
	groupname << "/reflevel=" << reflevel << "/";
	group = file.openGroup(groupname.str());

	// open Ye
	dataset = file.openDataSet(groupname.str()+"Ye");
	space = dataset.getSpace();
	comptype = dataset.getCompType();

	// get the dataset dimensions
	const int dataset_rank = 3;
	assert(space.getSimpleExtentNdims()==dataset_rank);                 // 3D array
	hsize_t space_dims[dataset_rank];
	space.getSimpleExtentDims(space_dims);
	nx = space_dims[0];
	ny = space_dims[1];
	nz = space_dims[2];
	int dataset_nzones = nx * ny * nz;

	//======================//
	// read the coordinates //
	//======================//
	// read in the deltas
	double delta[3];
	attr = group.openAttribute("delta");
	attr.read(H5::PredType::IEEE_F64LE,delta);
	for(int i=0; i<3; i++) delta[i] *= convert_length;
	dx = delta[0];
	dy = delta[1];
	dz = delta[2];
	assert(dx>0);
	assert(dy>0);
	assert(dz>0);

	// read in the extent.
	double extent[6];
	attr = group.openAttribute("extent");
	attr.read(H5::PredType::IEEE_F64LE,extent);
	for(int i=0; i<6; i++) extent[i] *= convert_length;

	// set the grid structure variables. half-cell offset since data is vertex-centered.
	x0 = extent[0] - dx/2.0;
	y0 = extent[2] - dy/2.0;
	z0 = extent[4] - dz/2.0;
	x1 = x0 + dx;
	y1 = y0 + dy;
	z1 = z0 + dz;
	xmax = extent[1] + dx/2.0;
	ymax = extent[3] + dy/2.0;
	zmax = extent[5] + dz/2.0;
	assert(xmax>x0);
	assert(ymax>y0);
	assert(zmax>z0);

	// check that the number of data points is consistent with the range and delta
	double ntmp;
	ntmp = (xmax - x0)/dx;
	assert(abs(nx-ntmp) < tiny);
	ntmp = (ymax - y0)/dy;
	assert(abs(ny-ntmp) < tiny);
	ntmp = (zmax - z0)/dz;
	assert(abs(nz-ntmp) < tiny);

	// modify grid in event of symmetry
	reflect_x = lua->scalar<int>("reflect_x");
	reflect_y = lua->scalar<int>("reflect_y");
	reflect_z = lua->scalar<int>("reflect_z");
	if(reflect_x){
		assert(xmax>0);
		x0 = 0.0;
		x1 = fmod(xmax,dx);
		if(abs(x1-x0)/dx < tiny) x1 += dx; // don't let the leftmost grid cell be stupidly tiny.
		assert(fmod(xmax-x1,dx) < tiny || fmod(xmax-x1,dx)-dx < tiny);
		nx = (int)((xmax-x1)/dx + 0.5) +1; // 0.5 to deal with finite precision.
	}
	if(reflect_y){
		assert(ymax>0);
		y0 = 0.0;
		y1 = fmod(ymax,dy);
		if(abs(y1-y0)/dy < tiny) y1 += dy; // don't let the leftmost grid cell be stupidly tiny.
		assert(fmod(ymax-y1,dy) < tiny || fmod(ymax-y1,dy)-dy < tiny);
		ny = (int)((ymax-y1)/dy + 0.5) +1; // 0.5 to deal with finite precision.
	}
	if(reflect_z){
		assert(zmax>0);
		z0 = 0.0;
		z1 = fmod(zmax,dz);
		if(abs(z1-z0)/dz < tiny) z1 += dz; // don't let the leftmost grid cell be stupidly tiny.
		assert(fmod(zmax-z1,dz) < tiny || fmod(zmax-z1,dz)-dz < tiny);
		nz = (int)((zmax-z1)/dz + 0.5) +1; // 0.5 to deal with finite precision.
	}

	// set up the zone structure
	int nzones = nx * ny * nz;
	z.resize(nzones,zone(3));

	// print out grid structure
	if(rank0){
		cout << "#   Using refinement level " << reflevel << endl;
		cout << "#   Minima           : {" << x0 << "," << y0 << "," << z0 <<"} cm" << endl;
		cout << "#   Next             : {" << x1 << "," << y1 << "," << z1 <<"} cm" << endl;
		cout << "#   Maxima           : {" << xmax << "," << ymax << "," << zmax <<"} cm" << endl;
		cout << "#   Deltas           : {" << dx << "," << dy << "," << dz <<"} cm" << endl;
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
	double Nmass=0.0, Rmass=0.0;
	double avgT=0, KE=0, totalVolume=0;
    #pragma omp parallel for collapse(3)
	for(unsigned i=0; i<space_dims[0]; i++)
		for(unsigned j=0; j<space_dims[1]; j++)
			for(unsigned k=0; k<space_dims[2]; k++){

				// get location of cell-centered cell center, i.e. vertex-centered vertex location
				vector<double> x(3,0.0);
				x[0] = extent[0] + i*dx;
				x[1] = extent[2] + j*dy;
				x[2] = extent[4] + k*dz;

				// get dataset and zone indices
				const int dataset_ind = k + space_dims[2]*j + space_dims[1]*space_dims[2]*i;
				const int z_ind = zone_index(x);
				assert(dataset_ind < (int)(space_dims[0] * space_dims[1] * space_dims[2]));
				assert(z_ind < (int)z.size());

				// if the grid cell is in our post-reflection-modified domain, add it to the zone list.
				if(z_ind >= 0){
					z[z_ind].rho  =  rho[dataset_ind];
					z[z_ind].T    = temp[dataset_ind];
					z[z_ind].Ye   =   Ye[dataset_ind];
					z[z_ind].v[0] = velx[dataset_ind];
					z[z_ind].v[1] = vely[dataset_ind];
					z[z_ind].v[2] = velz[dataset_ind];

					assert(z[z_ind].rho >= 0.0);
					assert(z[z_ind].T   >= 0.0);
					assert(z[z_ind].Ye >= 0.0);
					assert(z[z_ind].Ye <= 1.0);
					if(zone_speed2(z_ind) > pc::c*pc::c) cout << z_ind << " " << zone_speed2(z_ind)/pc::c/pc::c << endl;
					//assert(zone_speed2(z_ind) < pc::c*pc::c);
				}
			}

	// get global quantities
	#pragma omp parallel for reduction(+:Nmass,Rmass,avgT,KE,totalVolume)
	for(unsigned z_ind=0; z_ind<z.size(); z_ind++){
		Nmass += z[z_ind].rho * zone_lab_volume(z_ind);
		Rmass += z[z_ind].rho * zone_comoving_volume(z_ind);

		KE += 0.5 * z[z_ind].rho * zone_lab_volume(z_ind) * zone_speed2(z_ind);
		totalVolume += zone_lab_volume(z_ind);
		avgT += z[z_ind].T * zone_lab_volume(z_ind);
	}

	if(rank0){
		cout << "#   Newtonian    mass: " << Nmass * 2.0 << endl;
		cout << "#   Relativistic mass: " << Rmass * 2.0 << endl;
		cout << "#   Newtonian KE     : " << KE * 2.0 << endl;
		cout << "#   <T>              : " << avgT/totalVolume*pc::k_MeV << endl;
	}
	file.close();
}

void grid_3D_cart::read_SpEC_file(Lua* lua)
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
	infile >> x0;
	infile >> xmax;
	infile >> y0;
	infile >> ymax;
	infile >> z0;
	infile >> zmax;
	dx = (xmax-x0)/(double)nx; if(reflect_x) dx*=2;
	dy = (ymax-y0)/(double)ny; if(reflect_y) dy*=2;
	dz = (zmax-z0)/(double)nz; if(reflect_z) dz*=2;
	x1 = x0+dx;
	y1 = y0+dy;
	z1 = z0+dz;

	// First loop - set indices and read zone values in file
	// loop order is the file order
	// index order matches that in get_zone
	z.resize(n_zones,zone(dimensionality()));
	int ind = 0;
	bool rx,ry,rz;
	for (int k=0;k<nz;k++)
		for (int j=0;j<ny;j++)
			for (int i=0;i<nx;i++)
			{
				// set current index
				ind = zone_index(i,j,k);

				// create reverse map to x,y,z indices
				rx=false; ry=false; rz=false;
				if(reflect_x && i<nx/2) rx = true;
				if(reflect_y && j<ny/2) ry = true;
				if(reflect_z && k<nz/2) rz = true;

				// read in values if not in reflected zone
				if(!rx && !ry && !rz)
				{
					infile >> z[ind].rho;
					infile >> z[ind].T;
					infile >> z[ind].Ye;
					infile >> z[ind].v[0];
					infile >> z[ind].v[1];
					infile >> z[ind].v[2];
				}
				else{ //poison values
					z[ind].rho   = NaN;
					z[ind].T     = NaN;
					z[ind].Ye    = NaN;
					z[ind].v[0]  = NaN;
					z[ind].v[1]  = NaN;
					z[ind].v[2]  = NaN;
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
				ind = zone_index(i,j,k);

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
					origin_ind = zone_index(origin_i,origin_j,origin_k);
					z[ind].rho   = z[origin_ind].rho;
					z[ind].T = z[origin_ind].T;
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
int grid_3D_cart::zone_index(const vector<double>& x) const
{
	assert(x.size()==3);

	// check for off grid
	if (x[0]<x0 || x[0]>xmax) return -2;
	if (x[1]<y0 || x[1]>ymax) return -2;
	if (x[2]<z0 || x[2]>zmax) return -2;

	// get directional indices
	int i = (x[0]<=x1 ? 0 : floor((x[0]-x1)/dx)+1 );
	int j = (x[1]<=y1 ? 0 : floor((x[1]-y1)/dy)+1 );
	int k = (x[2]<=z1 ? 0 : floor((x[2]-z1)/dz)+1 );
	assert(i >= 0);
	assert(j >= 0);
	assert(k >= 0);
	assert(i < nx);
	assert(j < ny);
	assert(k < nz);

	int z_ind = zone_index(i,j,k);
	assert(z_ind >= 0);
	assert(z_ind < (int)z.size());
	return z_ind;
}


//-------------------------------------------------
// get the zone index from the directional indices
//-------------------------------------------------
int grid_3D_cart::zone_index(const int i, const int j, const int k) const{
	assert(i >= 0);
	assert(j >= 0);
	assert(k >= 0);
	assert(i < nx);
	assert(j < ny);
	assert(k < nz);
	const int z_ind = i*ny*nz + j*nz + k;
	assert(z_ind < (int)z.size());
	return z_ind;
}

//-------------------------------------------
// get directional indices from zone index
//-------------------------------------------
void grid_3D_cart::zone_directional_indices(const int z_ind, vector<int>& dir_ind) const
{
	assert(z_ind >= 0);
	assert(z_ind < (int)z.size());

	dir_ind.resize(dimensionality());
	dir_ind[0] =  z_ind / (ny*nz);
	dir_ind[1] = (z_ind % (ny*nz)) / nz;
	dir_ind[2] =  z_ind % nz;

	assert(dir_ind[0] >= 0);
	assert(dir_ind[0] < nx);
	assert(dir_ind[1] >= 0);
	assert(dir_ind[1] < ny);
	assert(dir_ind[2] >= 0);
	assert(dir_ind[2] < nz);
}


//------------------------------------------------------------
// return volume of zone (precomputed)
//------------------------------------------------------------
double grid_3D_cart::zone_lab_volume(const int z_ind) const
{
	assert(z_ind >= 0);
	assert(z_ind < (int)z.size());
	vector<int> dir_ind;
	zone_directional_indices(z_ind,dir_ind);
	vector<double> delta;
	get_deltas(z_ind,delta);
	return delta[0] * delta[1] * delta[2];
}


//------------------------------------------------------------
// sample a random position within the cubical cell
//------------------------------------------------------------
void grid_3D_cart::cartesian_sample_in_zone
(const int z_ind, const vector<double>& rand, vector<double>& x) const
{
	assert(z_ind >= 0);
	assert(z_ind < (int)z.size());
	x.resize(3);

	// zone directional indices
	vector<int> dir_ind;
	zone_directional_indices(z_ind,dir_ind);
	assert(dir_ind.size()==dimensionality());

	// zone deltas in each of three directions
	vector<double> delta;
	get_deltas(z_ind,delta);

	// set the random location
	x[0] = zone_left_boundary(0,dir_ind[0]) + delta[0]*rand[0];
	x[1] = zone_left_boundary(1,dir_ind[1]) + delta[1]*rand[1];
	x[2] = zone_left_boundary(2,dir_ind[2]) + delta[2]*rand[2];

	// make sure particle is in bounds
	assert(x[0] > x0 - tiny*dx);
	assert(x[1] > y0 - tiny*dy);
	assert(x[2] > z0 - tiny*dz);
	assert(x[0] < xmax + tiny*dx);
	assert(x[1] < ymax + tiny*dy);
	assert(x[2] < zmax + tiny*dz);

	// return particles just outside cell boundaries to within cell boundaries
	for(int i=0; i<3; i++){
		x[i] = min(x[i], zone_right_boundary(i,dir_ind[i]));
		x[i] = max(x[i],  zone_left_boundary(i,dir_ind[i]));
	}
}


//------------------------------------------------------------
// return length of zone
//------------------------------------------------------------
double  grid_3D_cart::zone_min_length(const int z_ind) const
{
	assert(z_ind >= 0);
	assert(z_ind < (int)z.size());

	vector<double> delta;
	get_deltas(z_ind,delta);

	double min_ds = min(delta[0], min(delta[1],delta[2]) );
	return min_ds;
}



//------------------------------------------------------------
// get the velocity vector 
//------------------------------------------------------------
void grid_3D_cart::cartesian_velocity_vector(const vector<double>& x, vector<double>& v, int z_ind) const
{
	assert(x.size()==3);
	v.resize(3);
	if(z_ind<0) z_ind = zone_index(x);
	assert(z_ind >= 0);
	assert(z_ind < (int)z.size());

	// may want to interpolate here?
	v[0] = z[z_ind].v[0];
	v[1] = z[z_ind].v[1];
	v[2] = z[z_ind].v[2];

	assert(v[0]*v[0] + v[1]*v[1] + v[2]*v[2] <= pc::c*pc::c);
}

//------------------------------------------------------------
// cell-centered coordinates of zone i
//------------------------------------------------------------
void grid_3D_cart::zone_coordinates(const int z_ind, vector<double>& r) const
{
	assert(z_ind >= 0);
	assert(z_ind < (int)z.size());
	r.resize(dimensionality());
	vector<int> dir_ind;
	zone_directional_indices(z_ind,dir_ind);
	r[0] = dir_ind[0]==0 ? 0.5*(x0+x1) : x1 + ((double)(dir_ind[0]-1) + 0.5) * dx;
	r[1] = dir_ind[1]==0 ? 0.5*(y0+y1) : y1 + ((double)(dir_ind[1]-1) + 0.5) * dy;
	r[2] = dir_ind[2]==0 ? 0.5*(z0+z1) : z1 + ((double)(dir_ind[2]-1) + 0.5) * dz;
}


//------------------------------------------------------------
// Write the grid information out to a file
//------------------------------------------------------------
void grid_3D_cart::write_rays(const int iw) const
{
	int i,j,k;
	int z_ind;
	vector<double> r;
	string filename = "";
	ofstream outf;

	// XY-slice
	filename = transport::filename("slice_xy",iw,".dat");
	outf.open(filename.c_str());
	write_header(outf);
	k = nz/2;
	for(i=0; i<nx; i++) for(j=0; j<ny; j++){
		if(j==0) outf << endl;
		z_ind = zone_index(i,j,k);
		zone_coordinates(z_ind,r);
		write_line(outf,z_ind);
	}
	outf.close();

	// XZ-slice
	filename = transport::filename("slice_xz",iw,".dat");
	outf.open(filename.c_str());
	write_header(outf);
	j = ny/2;
	for(i=0; i<nx; i++) for(k=0; k<nz; k++){
		if(k==0) outf << endl;
		z_ind = zone_index(i,j,k);
		zone_coordinates(z_ind,r);
		write_line(outf,z_ind);
	}
	outf.close();

	// YZ-slice
	filename = transport::filename("slice_yz",iw,".dat");
	outf.open(filename.c_str());
	write_header(outf);
	i = nx/2;
	for(j=0; j<ny; j++) for(k=0; k<nz; k++){
		if(k==0) outf << endl;
		z_ind = zone_index(i,j,k);
		zone_coordinates(z_ind,r);
		write_line(outf,z_ind);
	}
	outf.close();

	// X-direction
	filename = transport::filename("ray_x",iw,".dat");
	outf.open(filename.c_str());
	write_header(outf);
	j = ny/2;
	k = nz/2;
	for (i=0;i<nx;i++){
		z_ind = zone_index(i,j,k);
		zone_coordinates(z_ind,r);
		write_line(outf,z_ind);
	}
	outf.close();

	// Y-direction
	filename = transport::filename("ray_y",iw,".dat");
	outf.open(filename.c_str());
	write_header(outf);
	i = nx/2;
	k = nz/2;
	for (j=0; j<ny; j++){
		z_ind = zone_index(i,j,k);
		zone_coordinates(z_ind,r);
		write_line(outf,z_ind);
	}
	outf.close();

	// Z-direction
	filename = transport::filename("ray_z",iw,".dat");
	outf.open(filename.c_str());
	write_header(outf);
	i = nx/2;
	j = ny/2;
	for (k=0; k<nz; k++)
	{
		z_ind = zone_index(i,j,k);
		zone_coordinates(z_ind,r);
		write_line(outf,z_ind);
	}
	outf.close();
}


//------------------------------------------------------------
// Reflect off the outer boundary
//------------------------------------------------------------
void grid_3D_cart::reflect_outer(particle *p) const{
	// assumes particle is placed OUTSIDE of the zones
	int z_ind = zone_index(p->x);
	vector<double> delta;
	get_deltas(z_ind,delta);

	// invert the radial component of the velocity, put the particle just inside the boundary
	if(p->x[0] < x0){
		assert(p->D[0]<0);
		p->D[0] = -p->D[0];
		p->x[0] = x0 + tiny*delta[0];
	}
	if(p->x[1] < y0){
		assert(p->D[1]<0);
		p->D[1] = -p->D[1];
		p->x[2] = y0 + tiny*delta[1];
	}
	if(p->x[2] < z0){
		assert(p->D[0]<0);
		p->D[0] = -p->D[0];
		p->x[1] = z0 + tiny*delta[2];
	}
	if(p->x[0] > xmax){
		assert(p->D[0]>0);
		p->D[0] = -p->D[0];
		p->x[0] = xmax - tiny*delta[0];
	}
	if(p->x[1] > ymax){
		assert(p->D[1]>0);
		p->D[1] = -p->D[1];
		p->x[1] = ymax - tiny*delta[1];
	}
	if(p->x[2] > zmax){
		assert(p->D[2]>0);
		p->D[2] = -p->D[2];
		p->x[2] = zmax - tiny*delta[2];
	}

	// double check that the particle is in the boundary
	assert(p->x[0]>x0 && p->x[0]<xmax);
	assert(p->x[1]>y0 && p->x[1]<ymax);
	assert(p->x[2]>z0 && p->x[2]<zmax);
}


//------------------------------------------------------------
// Find distance to outer boundary
//------------------------------------------------------------
double grid_3D_cart::lab_dist_to_boundary(const particle *p) const{
	bool inside = (p->x[0] >= x0) && (p->x[0] <= xmax) &&
			(p->x[1] >= y0) && (p->x[1] <= ymax) &&
			(p->x[2] >= z0) && (p->x[2] <= zmax);

	// case: particle is inside the boundaries
	if(inside){
		double dist_x = (p->D[0]>0 ? xmax-p->x[0] : p->x[0]-x0);
		double dist_y = (p->D[1]>0 ? ymax-p->x[1] : p->x[1]-y0);
		double dist_z = (p->D[2]>0 ? zmax-p->x[2] : p->x[2]-z0);
		dist_x /= fabs(p->D[0]); assert(dist_x>=0);
		dist_y /= fabs(p->D[1]); assert(dist_y>=0);
		dist_z /= fabs(p->D[2]); assert(dist_z>=0);
		double dist = min(min(dist_x,dist_y),dist_z);
		assert(dist <= sqrt((xmax-x0)*(xmax-x0) + (ymax-y0)*(ymax-y0) + (zmax-z0)*(zmax-z0)));
		return dist;
	}

	// case: particle is outside the boundaries.
	else return numeric_limits<double>::infinity();
}


double grid_3D_cart::zone_radius(const int z_ind) const{
	assert(z_ind >= 0);
	assert(z_ind < (int)z.size());
	vector<double> r;
	zone_coordinates(z_ind,r);
	return sqrt(transport::dot(r,r));
}

//-----------------------------
// Dimensions of the grid
//-----------------------------
void grid_3D_cart::dims(vector<hsize_t>& dims) const{
	dims.resize(dimensionality());
	dims[0] = nx;
	dims[1] = ny;
	dims[2] = nz;
}

//----------------------------------------------------
// Write the coordinates of the grid points to the hdf5 file
//----------------------------------------------------
void grid_3D_cart::write_hdf5_coordinates(H5::H5File file) const
{
	// useful quantities
	H5::DataSet dataset;
	H5::DataSpace dataspace;
	vector<float> tmp;

	// get dimensions
	vector<hsize_t> coord_dims;
	dims(coord_dims);
	assert(coord_dims.size()==dimensionality());
	for(unsigned i=0; i<coord_dims.size(); i++) coord_dims[i]++; //make room for min value

	// write x coordinates
	dataspace = H5::DataSpace(1,&coord_dims[0]);
	dataset = file.createDataSet("grid_x(cm)",H5::PredType::IEEE_F32LE,dataspace);
	tmp.resize(coord_dims[0]);
	tmp[0] = x0;
	tmp[1] = x1;
	if(nx>1) for(int i=2; i<nx+1; i++) tmp[i] = x1 + (i-1)*dx;
	dataset.write(&tmp[0],H5::PredType::IEEE_F32LE);
	dataset.close();

	// write y coordinates
	dataspace = H5::DataSpace(1,&coord_dims[1]);
	dataset = file.createDataSet("grid_y(cm)",H5::PredType::IEEE_F32LE,dataspace);
	tmp.resize(coord_dims[1]);
	tmp[0] = y0;
	tmp[1] = y1;
	if(ny>1) for(int i=2; i<ny+1; i++) tmp[i] = y1 + (i-1)*dy;
	dataset.write(&tmp[0],H5::PredType::IEEE_F32LE);
	dataset.close();

	// write z coordinates
	dataspace = H5::DataSpace(1,&coord_dims[2]);
	dataset = file.createDataSet("grid_z(cm)",H5::PredType::IEEE_F32LE,dataspace);
	tmp.resize(coord_dims[2]);
	tmp[0] = z0;
	tmp[1] = z1;
	if(nz>1) for(int i=2; i<nz+1; i++) tmp[i] = z1 + (i-1)*dz;
	dataset.write(&tmp[0],H5::PredType::IEEE_F32LE);
	dataset.close();
}


void grid_3D_cart::get_deltas(const int z_ind, vector<double>& delta) const
{
	assert(z_ind < (int)z.size());
	delta.resize(dimensionality());

	// get directional indices
	vector<int> dir_ind;
	zone_directional_indices(z_ind,dir_ind);

	delta[0] = (dir_ind[0]==0 ? x1-x0 : dx);
	delta[1] = (dir_ind[1]==0 ? y1-y0 : dy);
	delta[2] = (dir_ind[2]==0 ? z1-z0 : dz);
	for(int i=0; i<3; i++) assert(delta[i]>0);
}


double grid_3D_cart::zone_left_boundary(const unsigned dir, const unsigned dir_ind) const{
	assert(dir>=0);
	assert(dir<3);
	assert(dir_ind>=0);

	double boundary = 0;
	if(dir==0)      boundary = ( dir_ind==0 ? x0 : x1 + (double)(dir_ind-1) * dx );
	else if(dir==1) boundary = ( dir_ind==0 ? y0 : y1 + (double)(dir_ind-1) * dy );
	else            boundary = ( dir_ind==0 ? z0 : z1 + (double)(dir_ind-1) * dz );

	return boundary;
}
double grid_3D_cart::zone_right_boundary(const unsigned dir, const unsigned dir_ind) const{
	assert(dir>=0);
	assert(dir<3);
	assert(dir_ind>=0);

	if(dir==0)      return dir_ind==0 ? x1 : x1 + (double)dir_ind * dx;
	else if(dir==1) return dir_ind==0 ? y1 : y1 + (double)dir_ind * dy;
	else            return dir_ind==0 ? z1 : z1 + (double)dir_ind * dz;

}
