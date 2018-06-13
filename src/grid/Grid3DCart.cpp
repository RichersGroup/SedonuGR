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
	PRINT_ASSERT(NDIMS,==,3);
	grid_type = "Grid3DCart";
	for(int i=0; i<3; i++) reflect[i] = 0;
	rotate_hemisphere[0] = 0;
	rotate_hemisphere[1] = 0;
	rotate_quadrant = 0;
	tetrad_rotation = cartesian;
}

void Grid3DCart::write_child_zones(H5::H5File file){
	if(DO_GR){
		betaup.write_HDF5(file,"shiftup");
		g3.write_HDF5(file,"threemetric");
		sqrtdetg3.write_HDF5(file,"sqrtdetg3");
	}
	v.write_HDF5(file,"threevelocity(cm|s)");
}

//------------------------------------------------------------
// Read in a cartesian model file
//------------------------------------------------------------
void Grid3DCart::read_model_file(Lua* lua)
{
	std::string model_type = lua->scalar<std::string>("model_type");
	if(model_type == "THC") read_THC_file(lua);
	else{
		cout << "ERROR: model type unknown." << endl;
		exit(8);
	}

	if(DO_GR){
		#pragma omp parallel for
		for(unsigned z_ind=0; z_ind<lapse.size(); z_ind++){
			Metric g;
			g.alpha = lapse[z_ind];
			for(unsigned i=0; i<3; i++) g.betaup[i] = betaup[z_ind][i];
			g.gammalow.data = g3[z_ind];
			g.update();
			sqrtdetg3[z_ind] = sqrt(g.gammalow.det());
		}
	}
}



void Grid3DCart::read_THC_file(Lua* lua)
{
	double dx[3], x0[3], x1[3], xmax[3];
	int nx[3];

	// verbocity
	int MPI_myID;
	MPI_Comm_rank( MPI_COMM_WORLD, &MPI_myID );
	const int rank0 = (MPI_myID == 0);

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
	int nzones = 1;
	for(int a=0; a<3; a++){
		vector<double> top(nx[a]), mid(nx[a]);
		top[0] = x1[a];
		mid[0] = 0.5 * (x0[a] + x1[a]);
		for(int i=1; i<nx[a]; i++){
			top[i] = x1[a] + (double)i*dx[a];
			mid[i] = 0.5 * (top[i] + top[i-1]);
		}
		xAxes[a] = Axis(x0[a], top, mid);
		nzones *= xAxes[a].size();
	}
	// set up the data structures
	v.set_axes(xAxes);
	rho.set_axes(xAxes);
	T.set_axes(xAxes);
	Ye.set_axes(xAxes);
	H_vis.set_axes(xAxes);
	if(DO_GR){
		betaup.set_axes(xAxes);
		g3.set_axes(xAxes);
		lapse.set_axes(xAxes);
		sqrtdetg3.set_axes(xAxes);
	}

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
	vector<double>   tmp_Ye(dataset_nzones,0.0);
	vector<double>  tmp_rho(dataset_nzones,0.0);
	vector<double> tmp_T(dataset_nzones,0.0);
	vector<double> tmp_velx(dataset_nzones,0.0);
	vector<double> tmp_vely(dataset_nzones,0.0);
	vector<double> tmp_velz(dataset_nzones,0.0);
	vector<double> tmp_betax(dataset_nzones,0.0);
	vector<double> tmp_betay(dataset_nzones,0.0);
	vector<double> tmp_betaz(dataset_nzones,0.0);
	vector<double> tmp_lapse(dataset_nzones,0.0);
	vector<double> tmp_gxx(dataset_nzones,0.0);
	vector<double> tmp_gxy(dataset_nzones,0.0);
	vector<double> tmp_gxz(dataset_nzones,0.0);
	vector<double> tmp_gyy(dataset_nzones,0.0);
	vector<double> tmp_gyz(dataset_nzones,0.0);
	vector<double> tmp_gzz(dataset_nzones,0.0);
	dataset = file.openDataSet(groupname.str()+"Ye");
	dataset.read(&tmp_Ye[0],H5::PredType::IEEE_F64LE);
	dataset = file.openDataSet(groupname.str()+"rho");
	dataset.read(&(tmp_rho[0]),H5::PredType::IEEE_F64LE);
	dataset = file.openDataSet(groupname.str()+"temp");
	dataset.read(&(tmp_T[0]),H5::PredType::IEEE_F64LE);
	dataset = file.openDataSet(groupname.str()+"velx");
	dataset.read(&(tmp_velx[0]),H5::PredType::IEEE_F64LE);
	dataset = file.openDataSet(groupname.str()+"vely");
	dataset.read(&(tmp_vely[0]),H5::PredType::IEEE_F64LE);
	dataset = file.openDataSet(groupname.str()+"velz");
	dataset.read(&(tmp_velz[0]),H5::PredType::IEEE_F64LE);
	if(DO_GR){
		dataset = file.openDataSet(groupname.str()+"betax");
		dataset.read(&(tmp_betax[0]),H5::PredType::IEEE_F64LE);
		dataset = file.openDataSet(groupname.str()+"betay");
		dataset.read(&(tmp_betay[0]),H5::PredType::IEEE_F64LE);
		dataset = file.openDataSet(groupname.str()+"betaz");
		dataset.read(&(tmp_betaz[0]),H5::PredType::IEEE_F64LE);
		dataset = file.openDataSet(groupname.str()+"lapse");
		dataset.read(&(tmp_lapse[0]),H5::PredType::IEEE_F64LE);
		dataset = file.openDataSet(groupname.str()+"gxx");
		dataset.read(&(tmp_gxx[0]),H5::PredType::IEEE_F64LE);
		dataset = file.openDataSet(groupname.str()+"gxy");
		dataset.read(&(tmp_gxy[0]),H5::PredType::IEEE_F64LE);
		dataset = file.openDataSet(groupname.str()+"gxz");
		dataset.read(&(tmp_gxz[0]),H5::PredType::IEEE_F64LE);
		dataset = file.openDataSet(groupname.str()+"gyy");
		dataset.read(&(tmp_gyy[0]),H5::PredType::IEEE_F64LE);
		dataset = file.openDataSet(groupname.str()+"gyz");
		dataset.read(&(tmp_gyz[0]),H5::PredType::IEEE_F64LE);
		dataset = file.openDataSet(groupname.str()+"gzz");
		dataset.read(&(tmp_gzz[0]),H5::PredType::IEEE_F64LE);
	}

	#pragma omp parallel for
	for(int i=0; i<dataset_nzones; i++){
		tmp_rho[i] *= convert_density;
		tmp_T[i] *= convert_temperature;
		tmp_velx[i] *= convert_velocity;
		tmp_vely[i] *= convert_velocity;
		tmp_velz[i] *= convert_velocity;
	}


	//===============//
	// fill the grid //
	//===============//
	#pragma omp parallel for
	for(int z_ind=0; z_ind<(int)rho.size(); z_ind++){

		// directional indices in Sedonu grid
		Tuple<unsigned,NDIMS> dir_ind = zone_directional_indices(z_ind);

		// directional indices in hdf5 data
		unsigned hdf5_dir_ind[3];
		for(int d=0; d<3; d++){
			hdf5_dir_ind[d] = offset[d] + dir_ind[d];
			PRINT_ASSERT(offset[d],>=,0);
		}

		// global hdf5 index
		const int dataset_ind = hdf5_dir_ind[2] + hdf5_dims[2]*hdf5_dir_ind[1] + hdf5_dims[1]*hdf5_dims[2]*hdf5_dir_ind[0];
		PRINT_ASSERT((int)dataset_ind,<,(int)(hdf5_dims[0] * hdf5_dims[1] * hdf5_dims[2]));
		PRINT_ASSERT((int)dataset_ind,>=,0);

		// fill the zone
		rho[z_ind]  =  tmp_rho[dataset_ind];
		T[z_ind]    =    tmp_T[dataset_ind];
		Ye[z_ind]   =   tmp_Ye[dataset_ind];
		v[z_ind][0] = tmp_velx[dataset_ind];
		v[z_ind][1] = tmp_vely[dataset_ind];
		v[z_ind][2] = tmp_velz[dataset_ind];
		if(DO_GR){
			lapse[z_ind]= tmp_lapse[dataset_ind];
			betaup[z_ind][0] = tmp_betax[dataset_ind];
			betaup[z_ind][1] = tmp_betay[dataset_ind];
			betaup[z_ind][2] = tmp_betaz[dataset_ind];
			g3[z_ind][ixx] = tmp_gxx[dataset_ind];
			g3[z_ind][ixy] = tmp_gxy[dataset_ind];
			g3[z_ind][ixz] = tmp_gxz[dataset_ind];
			g3[z_ind][iyy] = tmp_gyy[dataset_ind];
			g3[z_ind][iyz] = tmp_gyz[dataset_ind];
			g3[z_ind][izz] = tmp_gzz[dataset_ind];
		}

		PRINT_ASSERT(tmp_rho[z_ind],>=,0.0);
		PRINT_ASSERT(T[z_ind],>=,0.0);
		PRINT_ASSERT(tmp_Ye[z_ind],>=,0.0);
		PRINT_ASSERT(tmp_Ye[z_ind],<=,1.0);
	}

	file.close();
}

//------------------------------------------------------------
// Overly simple search to find zone
//------------------------------------------------------------
int Grid3DCart::zone_index(const Tuple<double,4>& x) const
{
	// check for off grid
	for(int i=0; i<3; i++){
		PRINT_ASSERT(x[i],>,xAxes[i].min   - xAxes[i].delta(1));
		PRINT_ASSERT(x[i],<,xAxes[i].max() + xAxes[i].delta(1));
		if (x[i]<xAxes[i].min || x[i]>=xAxes[i].max()) return -1;
	}

	// get directional indices
	int dir_ind[3] = {-1,-1,-1};
	for(int i=0; i<3; i++){
		dir_ind[i] = xAxes[i].bin(x[i]);
		PRINT_ASSERT(dir_ind[i],>=,0);
		PRINT_ASSERT(dir_ind[i],<,(int)xAxes[i].size());
	}

	int z_ind = zone_index(dir_ind[0],dir_ind[1],dir_ind[2]);
	PRINT_ASSERT(z_ind,>=,0);
	PRINT_ASSERT(z_ind,<,(int)rho.size());
	return z_ind;
}


//-------------------------------------------------
// get the zone index from the directional indices
//-------------------------------------------------
int Grid3DCart::zone_index(const int i, const int j, const int k) const{
	PRINT_ASSERT(i,>=,0);
	PRINT_ASSERT(j,>=,0);
	PRINT_ASSERT(k,>=,0);
	PRINT_ASSERT(i,<,(int)xAxes[0].size());
	PRINT_ASSERT(j,<,(int)xAxes[1].size());
	PRINT_ASSERT(k,<,(int)xAxes[2].size());
	const int z_ind = i*xAxes[1].size()*xAxes[2].size() + j*xAxes[2].size() + k;
	PRINT_ASSERT(z_ind,<,(int)rho.size());
	return z_ind;
}

//-------------------------------------------
// get directional indices from zone index
//-------------------------------------------
Tuple<unsigned,NDIMS> Grid3DCart::zone_directional_indices(const int z_ind) const
{
	PRINT_ASSERT(z_ind,>=,0);
	PRINT_ASSERT(z_ind,<,(int)rho.size());

	Tuple<unsigned,NDIMS> dir_ind;
	dir_ind[0] =  z_ind / (xAxes[1].size()*xAxes[2].size());
	dir_ind[1] = (z_ind % (xAxes[1].size()*xAxes[2].size())) / xAxes[2].size();
	dir_ind[2] =  z_ind % xAxes[2].size();

	for(int i=0; i<3; i++)
		PRINT_ASSERT(dir_ind[i],<,xAxes[i].size());
	return dir_ind;
}


//------------------------------------------------------------
// return volume of zone (precomputed)
//------------------------------------------------------------
double Grid3DCart::zone_lab_3volume(const int z_ind) const
{
	PRINT_ASSERT(z_ind,>=,0);
	PRINT_ASSERT(z_ind,<,(int)rho.size());
	double delta[3];
	get_deltas(z_ind,delta,3);
	double result = delta[0] * delta[1] * delta[2];
	if(DO_GR) result *= sqrtdetg3[z_ind];
	PRINT_ASSERT(result,>,0);
	return result;
}


//------------------------------------------------------------
// sample a random position within the cubical cell
//------------------------------------------------------------
Tuple<double,4> Grid3DCart::sample_in_zone(const int z_ind, ThreadRNG* rangen) const
{
	PRINT_ASSERT(z_ind,>=,0);
	PRINT_ASSERT(z_ind,<,(int)rho.size());

	double rand[3];
	rand[0] = rangen->uniform();
	rand[1] = rangen->uniform();
	rand[2] = rangen->uniform();

	// zone directional indices
	Tuple<unsigned,NDIMS> dir_ind = zone_directional_indices(z_ind);

	// zone deltas in each of three directions
	double delta[3];
	get_deltas(z_ind,delta,3);

	// set the random location
	Tuple<double,4> x;
	for(int i=0; i<3; i++){
		x[i] = zone_left_boundary(i,dir_ind[i]) + delta[i]*rand[i];

		// make sure the particle is in bounds
		PRINT_ASSERT(x[i],>,xAxes[i].min   - TINY*xAxes[i].delta(0));
		PRINT_ASSERT(x[0],<,xAxes[i].max() + TINY*xAxes[i].delta(0));

		// return particles just outside cell boundaries to within cell boundaries
		x[i] = min(x[i], zone_right_boundary(i,dir_ind[i]));
		x[i] = max(x[i],  zone_left_boundary(i,dir_ind[i]));
	}
	return x;
}


//------------------------------------------------------------
// return length of zone
//------------------------------------------------------------
double  Grid3DCart::zone_min_length(const int z_ind) const
{
	PRINT_ASSERT(z_ind,>=,0);
	PRINT_ASSERT(z_ind,<,(int)rho.size());

	double delta[3];
	get_deltas(z_ind,delta,3);

	double min_ds = min(delta[0], min(delta[1],delta[2]) );
	return min_ds;
}

// returning 0 causes the min distance to take over in propagate.cpp::which_event
double Grid3DCart::d_boundary(const EinsteinHelper *eh) const{

	// x direction
	double dlambda[3] = {INFINITY,INFINITY,INFINITY};
	for(unsigned d=0; d<3; d++){
		unsigned i = eh->dir_ind[d];
		PRINT_ASSERT(eh->xup[d],<=,zone_right_boundary(d,i));
		PRINT_ASSERT(eh->xup[d],>=, zone_left_boundary(d,i));
		if(eh->kup[d] < 0) dlambda[d] = ( zone_left_boundary(d,i) - eh->xup[d]) / eh->kup[d];
		if(eh->kup[d] > 0) dlambda[d] = (zone_right_boundary(d,i) - eh->xup[d]) / eh->kup[d];
		PRINT_ASSERT(dlambda[d],>=,0);
	}

	double ds_com = min( dlambda[0], min(dlambda[1], dlambda[2])) * eh->kup_tet[3];
	PRINT_ASSERT(ds_com,>=,0);
	return ds_com;
}

double Grid3DCart::d_randomwalk(const EinsteinHelper *eh) const{
	double R=INFINITY;
	double D = eh->scatopac / (3.*pc::c);

	for(unsigned i=0; i<3; i++){
		for(int sgn=1; sgn>0; sgn*=-1){
			// get a null test vector
			Tuple<double,4> ktest;
			for(size_t j=0; j<4; j++) ktest[j] = 0;
			ktest[i] = sgn;
			eh->g.normalize_null_changeupt(ktest);
			if(ktest[3]<0) for(unsigned i=0; i<4; i++) ktest[i] *= -1;

			// get the time component of the tetrad test vector
			double kup_tet_t = -eh->g.dot<4>(ktest,eh->u);
			PRINT_ASSERT(kup_tet_t,>,0);

			// get the min distance from the boundary in direction i. Negative if moving left
			double dxlab=0;
			if(sgn>0) dxlab = xAxes[i].top[eh->dir_ind[i]] - eh->xup[i];
			if(sgn<0) dxlab = xAxes[i].bottom(eh->dir_ind[i]) - eh->xup[i];

			R = min(R, sim->R_randomwalk(ktest[i]/kup_tet_t, ktest[3]/kup_tet_t, eh->u[i], dxlab, D));
		}
	}

	PRINT_ASSERT(R,>=,0);
	PRINT_ASSERT(R,<,INFINITY);
	return R;
}
//------------------------------------------------------------
// get the velocity vector 
//------------------------------------------------------------
Tuple<double,3> Grid3DCart::interpolate_fluid_velocity(const EinsteinHelper& eh) const
{
	// may want to interpolate here?
	return v.interpolate(eh.icube_vol);
}

//------------------------------------------------------------
// cell-centered coordinates of zone i
//------------------------------------------------------------
Tuple<double,NDIMS> Grid3DCart::zone_coordinates(const int z_ind) const
{
	PRINT_ASSERT(z_ind,>=,0);
	PRINT_ASSERT(z_ind,<,(int)rho.size());
	Tuple<unsigned,NDIMS> dir_ind = zone_directional_indices(z_ind);

	Tuple<double,NDIMS> r;
	for(int i=0; i<3; i++)
		r[i] = xAxes[i].mid[dir_ind[i]];
	return r;
}


double Grid3DCart::zone_radius(const int z_ind) const{
	PRINT_ASSERT(z_ind,>=,0);
	PRINT_ASSERT(z_ind,<,(int)rho.size());
	Tuple<double,NDIMS> r = zone_coordinates(z_ind);
	return radius(r);
}

//-----------------------------
// Dimensions of the grid
//-----------------------------
Tuple<hsize_t,NDIMS> Grid3DCart::dims() const{
	Tuple<hsize_t,NDIMS> dims;
	for(int i=0; i<3; i++) dims[i] = xAxes[i].size();
	return dims;
}

void Grid3DCart::get_deltas(const int z_ind, double delta[3], const int size) const
{
	PRINT_ASSERT(z_ind,<,(int)rho.size());
	PRINT_ASSERT(size,==,3);

	// get directional indices
	Tuple<unsigned,NDIMS> dir_ind = zone_directional_indices(z_ind);

	for(int i=0; i<3; i++){
		delta[i] =xAxes[i].delta(dir_ind[i]);
		PRINT_ASSERT(delta[i],>,0);
	}
}


double Grid3DCart::zone_left_boundary(const unsigned dir, const unsigned dir_ind) const{
	PRINT_ASSERT(dir,<,3);

	double boundary = xAxes[dir].bottom(dir_ind);
	PRINT_ASSERT(boundary,<=,xAxes[dir].max());
	PRINT_ASSERT(boundary,>=,xAxes[dir].min);
	return boundary;
}
double Grid3DCart::zone_right_boundary(const unsigned dir, const unsigned dir_ind) const{
	PRINT_ASSERT(dir,<,3);

	double boundary = xAxes[dir].top[dir_ind];
	PRINT_ASSERT(boundary,<=,xAxes[dir].max()*(1.0+TINY));
	PRINT_ASSERT(boundary,>=,xAxes[dir].min);
	return boundary;
}


//------------------------------------------------------------
// Reflect off revlecting boundary condition
//------------------------------------------------------------
void Grid3DCart::symmetry_boundaries(EinsteinHelper *eh) const{
	PRINT_ASSERT(eh->fate,==,moving);

	// initialize the arrays
	double kup[4], xup[4];
	for(int i=0; i<4; i++) xup[i] = eh->xup[i];
	for(int i=0; i<4; i++) kup[i] = eh->kup[i];


	// invert the radial component of the velocity, put the particle just inside the boundary
	for(int i=0; i<3; i++){
		if(reflect[i] && xup[i]<0){
			PRINT_ASSERT(xAxes[i].min,==,0);
			PRINT_ASSERT(-xup[i],<,xAxes[i].delta(1));
			// actual work
			kup[i] = -kup[i];
			xup[i] = -xup[i];
			// end actual work
			PRINT_ASSERT(xup[i],>=,xAxes[i].min);
			PRINT_ASSERT(xup[i],<=,xAxes[i].max());
		}
	}

	// rotating boundary conditions
	for(int i=0; i<2; i++){
		if(xup[i]<0 && (rotate_hemisphere[i] || rotate_quadrant)){
			PRINT_ASSERT(xAxes[i].min,==,0);
			PRINT_ASSERT(-xup[i],<,xAxes[i].delta(1));
			
			if(rotate_hemisphere[i]){
				for(int j=0; j<2; j++){
					kup[j] = -kup[j];
					xup[j] = -xup[j];
				}
			}
			else if(rotate_quadrant){
				int other = i==0 ? 1 : 0;
				if(xup[other]>=0){
					double tmp;
					tmp=kup[i];	kup[i]=kup[other]; kup[other]=-tmp;
					tmp=xup[i];	xup[i]=xup[other]; xup[other]=-tmp;
				}
				else for(int j=0; j<2; j++){
					kup[j] = -kup[j];
					xup[j] = -xup[j];
				}
			}

			// double check that the particle is in the boundary
			PRINT_ASSERT(xup[0],>=,xAxes[0].min);
			PRINT_ASSERT(xup[1],>=,xAxes[1].min);
			//PRINT_ASSERT(xup[0],<=,xAxes[0].max());
			//PRINT_ASSERT(xup[1],<=,xAxes[1].max());
		}
	}

	// assign the arrays
	for(unsigned i=0; i<4; i++){
		eh->xup[i] = xup[i];
		eh->kup[i] = kup[i];
	}
}

double Grid3DCart::zone_lorentz_factor(const int z_ind) const{
	Metric g;
	if(DO_GR) g.gammalow.data = g3[z_ind];
	else g.gammalow.data = NaN;

	Tuple<double,3> threevel;
	for(unsigned i=0; i<3; i++) threevel[i] = v[z_ind][i]/pc::c;

	double result = EinsteinHelper::lorentzFactor(&g,threevel);
	PRINT_ASSERT(result,>=,1);
	PRINT_ASSERT(result,<,INFINITY);
	return result;
}

Tuple<double,4> Grid3DCart::dk_dlambda(const EinsteinHelper& eh) const{
  double dg[4][4][4];
  for(unsigned i=0; i<4; i++) for(unsigned j=0; j<4; j++){ // no time derivatives
      dg[3][i][j] = 0;
    }

  Tuple<Tuple<double,6>,NDIMS> dg3_dx = g3.interpolate_slopes(eh.icube_vol);
  Tuple<double,NDIMS> da_dx = lapse.interpolate_slopes(eh.icube_vol);
  Tuple<Tuple<double,3>,NDIMS> dbetaup_dx = betaup.interpolate_slopes(eh.icube_vol);
  Tuple<Tuple<double,3>,NDIMS> dbetalow_dx;

  for(unsigned a=0; a<3; a++) 
    dbetalow_dx[a] = eh.g.gammalow.lower(dbetaup_dx[a]);
  
  #pragma omp simd
  for(unsigned a=0; a<3; a++){
    // xx parts
    dg[a][0][0]             = dg3_dx[a][ixx]; // [direction][element]
    dg[a][1][1]             = dg3_dx[a][iyy];
    dg[a][2][2]             = dg3_dx[a][izz];
    dg[a][0][1]=dg[a][1][0] = dg3_dx[a][ixy];
    dg[a][0][2]=dg[a][2][0] = dg3_dx[a][ixz];
    dg[a][1][2]=dg[a][2][1] = dg3_dx[a][iyz];

    // xt and tt parts
    dg[a][3][3] = -2.*eh.g.alpha*da_dx[a];
    for(unsigned i=0; i<3; i++){
      dg[a][3][3] += 2. * eh.g.betalow[i] * dbetaup_dx[a][i]; // [direction][element]
      dg[a][3][i] = dbetalow_dx[a][i];
      dg[a][i][3] = dbetalow_dx[a][i];
    }
  }

  // get the low-index Christoffel symbols
  Tuple<double,4> dk_dlambda_low = 0;
  #pragma omp simd collapse(3)
  for(unsigned a=0; a<4; a++)
    for(unsigned i=0; i<4; i++)
      for(unsigned j=0; j<4; j++)
	dk_dlambda_low[a] += (dg[i][a][j] - 0.5*dg[a][i][j]) * eh.kup[i] * eh.kup[j];

  Tuple<double,4> dk_dlambda = eh.g.raise(dk_dlambda_low);
  return dk_dlambda * -1.;
}
void Grid3DCart::interpolate_shift(EinsteinHelper* eh) const{
	if(DO_GR){
		Tuple<double,3> tmp = betaup.interpolate(eh->icube_vol);
		for(unsigned i=0; i<3; i++) eh->g.betaup[i] = tmp[i];
	}
}
void Grid3DCart::interpolate_3metric(EinsteinHelper* eh) const{
	if(DO_GR) eh->g.gammalow.data = g3.interpolate(eh->icube_vol);
}
void Grid3DCart::grid_coordinates(const Tuple<double,4>& xup, double coords[NDIMS]) const{
	coords[0] = xup[0];
	coords[1] = xup[1];
	coords[2] = xup[2];
}
