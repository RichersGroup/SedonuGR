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
#include "Transport.h"
#include "Grid1DSphere.h"
#include "global_options.h"
#include <cmath>

using namespace std;
namespace pc = physical_constants;

Grid1DSphere::Grid1DSphere(){
	PRINT_ASSERT(NDIMS,==,1);
	grid_type = "Grid1DSphere";
	reflect_outer = 0;
	tetrad_rotation = spherical;
}

//------------------------------------------------------------
// initialize the zone geometry from model file
//------------------------------------------------------------
void Grid1DSphere::read_model_file(Lua* lua)
{
	std::string model_type = lua->scalar<std::string>("model_type");
	if(model_type == "Nagakura") read_nagakura_model(lua);
	else if(model_type == "custom") read_custom_model(lua);
	else{
		cout << "ERROR: model type unknown." << endl;
		exit(8);
	}

	grid_type = "Grid1DSphere";
	reflect_outer = lua->scalar<int>("reflect_outer");

	vr.calculate_slopes(-INFINITY,INFINITY);
	X.calculate_slopes(0,INFINITY);
}

void Grid1DSphere::read_nagakura_model(Lua* lua){
	// verbocity
	int MPI_myID;
	MPI_Comm_rank( MPI_COMM_WORLD, &MPI_myID );
	const int rank0 = (MPI_myID == 0);
	vector<double> bintops, binmid;
	double trash, minval, tmp;

	// open the model files
	if(rank0) cout << "# Reading the model file..." << endl;
	string model_file = lua->scalar<string>("model_file");
	ifstream infile;
	infile.open(model_file.c_str());
	if(infile.fail()){
		if(rank0) cout << "Error: can't read the model file." << model_file << endl;
		exit(4);
	}


	// read in the radial grid
	string rgrid_filename = lua->scalar<string>("Grid1DSphere_Nagakura_rgrid_file");
	ifstream rgrid_file;
	rgrid_file.open(rgrid_filename.c_str());
	rgrid_file >> trash >> minval;
	bintops = vector<double>(0);
	binmid = vector<double>(0);
	while(rgrid_file >> trash >> tmp){
		if(bintops.size()>0) PRINT_ASSERT(tmp,>,bintops[bintops.size()-1]);
		else PRINT_ASSERT(tmp,>,minval);
		bintops.push_back(tmp);
		double last = bintops.size()==1 ? minval : bintops[bintops.size()-2];
		double midpoint = 0.5 * (bintops[bintops.size()-1] + last);
		binmid.push_back(midpoint);
	}
	rAxis = Axis(minval, bintops, binmid);
	vector<Axis> axes = {rAxis};
	vr.set_axes(axes);
	lapse.set_axes(axes);
	X.set_axes(axes);
	rho.set_axes(axes);
	T.set_axes(axes);
	Ye.set_axes(axes);
	H_vis.set_axes(axes);

	// write grid properties
	if(rank0)
	  cout << "#   nr=" << rAxis.size() << "\trmin=" << rAxis.min << "\trmax=" << rAxis.top[rAxis.size()-1] << endl;

	// read the fluid properties
	for(int z_ind=0; z_ind<rAxis.size(); z_ind++){
		double trash;

		// read the contents of a single line
		infile >> trash; // r
		infile >> trash; // theta
		infile >> rho[z_ind]; // g/ccm
		infile >> Ye[z_ind];
		infile >> T[z_ind]; // MeV
		infile >> vr[z_ind]; // cm/s
		infile >> trash; // 1/s
		infile >> trash; // 1/s

		// get rid of the rest of the line
		for(int k=9; k<=165; k++) infile >> trash;

		// convert units
		vr[z_ind] *= rAxis.mid[z_ind];
		T[z_ind] /= pc::k_MeV;

		// sanity checks
		PRINT_ASSERT(rho[z_ind],>=,0.0);
		PRINT_ASSERT(T[z_ind],>=,0.0);
		PRINT_ASSERT(Ye[z_ind],>=,0.0);
		PRINT_ASSERT(Ye[z_ind],<=,1.0);
	}
}

void Grid1DSphere::read_custom_model(Lua* lua){
	// verbocity
	int MPI_myID;
	MPI_Comm_rank( MPI_COMM_WORLD, &MPI_myID );
	const int rank0 = (MPI_myID == 0);
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
	PRINT_ASSERT(n_zones,>,0);
	vector<double> rtop(n_zones), rmid(n_zones);
	double rmin;

	// read zone properties
	vector<double> tmp_rho = vector<double>(n_zones,0);
	vector<double> tmp_T = vector<double>(n_zones,0);
	vector<double> tmp_Ye = vector<double>(n_zones,0);
	vector<double> tmp_H_vis = vector<double>(n_zones,0);
	vector<double> tmp_alpha = vector<double>(n_zones,0);
	vector<double> tmp_X = vector<double>(n_zones,0);
	vector<double> tmp_vr = vector<double>(n_zones,0);
	infile >> rmin;
	PRINT_ASSERT(rmin,>=,0);
	for(int z_ind=0; z_ind<n_zones; z_ind++)
	{
		double alpha, x;
		infile >> rtop[z_ind];
		infile >> tmp_rho[z_ind];
		infile >> tmp_T[z_ind];
		infile >> tmp_Ye[z_ind];
		tmp_H_vis[z_ind] = 0;
		infile >> tmp_vr[z_ind];
		infile >> tmp_alpha[z_ind];
		infile >> tmp_X[z_ind];

		double last = z_ind==0 ? rmin : rtop[z_ind-1];
		rmid[z_ind] = 0.5 * (rtop[z_ind] + last);
		PRINT_ASSERT(rtop[z_ind],>,(z_ind==0 ? rmin : rtop[z_ind-1]));
		PRINT_ASSERT(tmp_rho[z_ind],>=,0);
		PRINT_ASSERT(tmp_T[z_ind],>=,0);
		PRINT_ASSERT(tmp_Ye[z_ind],>=,0);
		PRINT_ASSERT(tmp_Ye[z_ind],<=,1.0);
		PRINT_ASSERT(tmp_alpha[z_ind],<=,1.0);
		PRINT_ASSERT(tmp_X[z_ind],>=,1.0);
	}
	rAxis = Axis(rmin, rtop, rmid);
	vector<Axis> axes = {rAxis};
	vr.set_axes(axes);
	lapse.set_axes(axes);
	X.set_axes(axes);
	rho.set_axes(axes);
	T.set_axes(axes);
	Ye.set_axes(axes);
	H_vis.set_axes(axes);

	for(unsigned z_ind=0; z_ind<vr.size(); z_ind++){
		vr[z_ind] = tmp_vr[z_ind];
		lapse[z_ind] = tmp_alpha[z_ind];
		X[z_ind] = tmp_X[z_ind];
		rho[z_ind] = tmp_rho[z_ind];
		T[z_ind] = tmp_T[z_ind];
		Ye[z_ind] = tmp_Ye[z_ind];
		H_vis[z_ind] = tmp_H_vis[z_ind];
	}

	// set the christoffel symbol coefficients
	vector<double> tmp_dadr = vector<double>(n_zones,0);
	vector<double> tmp_dXdr = vector<double>(n_zones,0);

	infile.close();
}


//------------------------------------------------------------
// Return the zone index containing the position x
//------------------------------------------------------------
int Grid1DSphere::zone_index(const double x[3]) const
{
	PRINT_ASSERT(rho.size(),>,0);
	double r = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
	PRINT_ASSERT(r,>=,0);

	// check if off the boundaries
	if(r < rAxis.min             ) return -1;
	if(r >= rAxis.top[rAxis.size()-1] ) return -1;

	// find in zone array using stl algorithm upper_bound and subtracting iterators
	int z_ind = rAxis.bin(r);
	PRINT_ASSERT(z_ind,>=,0);
	PRINT_ASSERT(z_ind,<,(int)rho.size());
	return z_ind;
}


//------------------------------------------------------------
// return volume of zone z_ind
//------------------------------------------------------------
double  Grid1DSphere::zone_lab_3volume(const int z_ind) const
{
	PRINT_ASSERT(z_ind,>=,0);
	PRINT_ASSERT(z_ind,<,(int)rho.size());
	double r0 = (z_ind==0 ? rAxis.min : rAxis.top[z_ind-1]);
	double vol = 4.0*pc::pi/3.0*( pow(rAxis.top[z_ind],3) - pow(r0,3) );
	if(DO_GR) vol *= X[z_ind];
	PRINT_ASSERT(vol,>=,0);
	return vol;
}

//------------------------------------------------------------
// return length of zone
//------------------------------------------------------------
double  Grid1DSphere::zone_min_length(const int z_ind) const
{
	PRINT_ASSERT(z_ind,>=,0);
	PRINT_ASSERT(z_ind,<,(int)rho.size());
	double r0 = (z_ind==0 ? rAxis.min : rAxis.top[z_ind-1]);
	double min_len = rAxis.top[z_ind] - r0;
	PRINT_ASSERT(min_len,>=,0);
	return min_len;
}


// ------------------------------------------------------------
// find the coordinates of the zone in geometrical coordinates
// ------------------------------------------------------------
void Grid1DSphere::zone_coordinates(const int z_ind, double r[1], const int rsize) const{
	PRINT_ASSERT(z_ind,>=,0);
	PRINT_ASSERT(z_ind,<,(int)rho.size());
	PRINT_ASSERT(rsize,==,(int)dimensionality());
	r[0] = 0.5*(rAxis.top[z_ind]+rAxis.bottom(z_ind));
	PRINT_ASSERT(r[0],>,0);
	PRINT_ASSERT(r[0],<,rAxis.top[rAxis.size()-1]);
}


//-------------------------------------------
// get directional indices from zone index
//-------------------------------------------
void Grid1DSphere::zone_directional_indices(const int z_ind, vector<unsigned>& dir_ind) const
{
	PRINT_ASSERT(z_ind,>=,0);
	PRINT_ASSERT(z_ind,<,(int)rho.size());
	PRINT_ASSERT(dir_ind.size(),==,(int)dimensionality());
	dir_ind[0] = z_ind;
}


//------------------------------------------------------------
// sample a random position within the spherical shell
//------------------------------------------------------------
void Grid1DSphere::sample_in_zone(const int z_ind, ThreadRNG* rangen, double x[3]) const
{
	PRINT_ASSERT(z_ind,>=,0);
	PRINT_ASSERT(z_ind,<,(int)rho.size());

	double rand[3];
	rand[0] = rangen->uniform();
	rand[1] = rangen->uniform();
	rand[2] = rangen->uniform();

	// inner and outer radii of shell
	double r0 = (z_ind==0 ? rAxis.min : rAxis.top[z_ind-1]);
	double r1 = rAxis.top[z_ind];

	// sample radial position in shell using a probability integral transform
	double radius = pow( rand[0]*(r1*r1*r1 - r0*r0*r0) + r0*r0*r0, 1./3.);
	PRINT_ASSERT(radius,>=,r0*(1.-TINY));
	PRINT_ASSERT(radius,<=,r1*(1.+TINY));
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
void Grid1DSphere::interpolate_fluid_velocity(const double x[3], double v[3], const unsigned dir_ind[NDIMS]) const
{
	// radius in zone
	double r = radius(x);

	// assuming radial velocity (may want to interpolate here)
	// (the other two components are ignored and mean nothing)
	double vr_interp = vr.interpolate(&r,dir_ind);
	v[0] = x[0]/r*vr_interp;
	v[1] = x[1]/r*vr_interp;
	v[2] = x[2]/r*vr_interp;

	// check for pathological case
	if (r == 0)
	{
		v[0] = 0;
		v[1] = 0;
		v[2] = 0;
	}

	PRINT_ASSERT(Metric::dot_Minkowski<3>(v,v),<=,pc::c*pc::c);
}



//------------------------------------------------------------
// Reflect off symmetry axis
//------------------------------------------------------------
void Grid1DSphere::symmetry_boundaries(EinsteinHelper *eh) const{
	// reflect from outer boundary
	double R = radius(eh->p.xup);
	if(reflect_outer && R>rAxis.top[rAxis.size()-1]){
		double r0 = (rAxis.size()>1 ? rAxis.top[rAxis.size()-2] : rAxis.min);
		double rmax = rAxis.top[rAxis.size()-1];
		double dr = rmax - r0;

		double kr;
		for(int i=0; i<3; i++) kr += eh->p.xup[i]/R * eh->p.kup[i];

		// invert the radial component of the velocity
		eh->p.kup[0] -= 2.*kr * eh->p.xup[0]/R;
		eh->p.kup[1] -= 2.*kr * eh->p.xup[1]/R;
		eh->p.kup[2] -= 2.*kr * eh->p.xup[2]/R;
		eh->g.normalize_null(eh->p.kup);

		// put the particle just inside the boundary
		double newR = rmax - TINY*dr;
		for(unsigned i=0; i<3; i++)	eh->p.xup[i] *= newR/R;

		// must be inside the boundary, or will get flagged as escaped
		PRINT_ASSERT(zone_index(eh->p.xup),>=,0);
	}
}

double Grid1DSphere::d_boundary(const EinsteinHelper *eh) const{
	double r = radius(eh->p.xup);
	PRINT_ASSERT(r,<=,rAxis.top[eh->z_ind]);
	PRINT_ASSERT(r,>=,rAxis.bottom(eh->z_ind));
	const double x = eh->p.xup[0];
	const double y = eh->p.xup[1];
	const double z = eh->p.xup[2];

	// get component of k in the radial direction
	double kr = eh->g.dot<4>(eh->e[0],eh->p.kup);

	double dlambda = INFINITY;
	if(kr>0) dlambda = (rAxis.top[eh->z_ind] - r   ) / kr;
	if(kr<0) dlambda = (r = rAxis.bottom(eh->z_ind)) / kr;

	return dlambda * eh->kup_tet[3];
}
double Grid1DSphere::d_randomwalk(const EinsteinHelper *eh) const{
	double R=INFINITY;
	double D = eh->scatopac / (3.*pc::c);

	double ktest[4] = {eh->p.xup[0], eh->p.xup[1], eh->p.xup[2], 0};
	const double r = radius(eh->p.xup);
	const double kr = r;
	const double ur = radius(eh->u);

	for(int sgn=1; sgn>0; sgn*=-1){
		// get a null test vector
		for(unsigned i=0; i<3; i++) ktest[i] *= sgn;
		eh->g.normalize_null(ktest);

		// get the time component of the tetrad test vector
		double kup_tet_t = -eh->g.dot<4>(ktest,eh->u);

		// get the min distance from the boundary in direction i. Negative if moving left
		double drlab=0;
		if(sgn>0) drlab = rAxis.top[eh->dir_ind[0]] - r;
		if(sgn<0) drlab = rAxis.bottom(eh->dir_ind[0]) - r;

		R = min(R, sim->R_randomwalk(kr/kup_tet_t, ur, drlab, D));
	}

	PRINT_ASSERT(R,>=,0);
	PRINT_ASSERT(R,<,INFINITY);
	return R;
}

double Grid1DSphere::zone_radius(const int z_ind) const{
	PRINT_ASSERT(z_ind,>=,0);
	PRINT_ASSERT(z_ind,<,(int)rho.size());
	return rAxis.top[z_ind];
}

//-----------------------------
// Dimensions of the grid
//-----------------------------
void Grid1DSphere::dims(hsize_t dims[1], const int size) const{
	PRINT_ASSERT(size,==,(int)dimensionality());
	dims[0] = rAxis.size();
}

//----------------------------------------------------
// Write the coordinates of the grid points to the hdf5 file
//----------------------------------------------------
void Grid1DSphere::write_hdf5_coordinates(H5::H5File file) const
{
	// useful quantities
	H5::DataSet dataset;
	H5::DataSpace dataspace;
	float tmp[rAxis.size()+1];

	// get dimensions
	hsize_t coord_dims[1];
	dims(coord_dims,1);
	for(unsigned i=0; i<1; i++) coord_dims[i]++; //make room for min value

	// write coordinates
	dataspace = H5::DataSpace(1,&coord_dims[0]);
	dataset = file.createDataSet("r(cm)",H5::PredType::IEEE_F32LE,dataspace);
	tmp[0] = rAxis.min;
	for(unsigned i=1; i<rAxis.size()+1; i++) tmp[i] = rAxis.top[i-1];
	dataset.write(&tmp[0],H5::PredType::IEEE_F32LE);
	dataset.close();
}
void Grid1DSphere::interpolate_3metric(const double xup[4], ThreeMetric* gammalow, const unsigned dir_ind[NDIMS]) const{
	const double r = radius(xup);
	const double Xloc = X.interpolate(&r,dir_ind);//1./sqrt(1.-1./r);//
	double tmp = (Xloc*Xloc-1.0) / (r*r);

	gammalow->data[ixx] = xup[0]*xup[0] * tmp;
	gammalow->data[iyy] = xup[1]*xup[1] * tmp;
	gammalow->data[izz] = xup[2]*xup[2] * tmp;
	gammalow->data[ixy] = xup[0]*xup[1] * tmp;
	gammalow->data[ixz] = xup[0]*xup[2] * tmp;
	gammalow->data[iyz] = xup[1]*xup[2] * tmp;

	gammalow->data[ixx] += 1.0;
	gammalow->data[iyy] += 1.0;
	gammalow->data[izz] += 1.0;
}

void Grid1DSphere::get_connection_coefficients(EinsteinHelper* eh) const{
	const double r = radius(eh->p.xup);
	const double alpha = lapse.interpolate(eh->p.xup,eh->dir_ind); //sqrt(1.-1./r); //
	const double Xloc = X.interpolate(eh->p.xup, eh->dir_ind); //1./alpha; //
	const double dadr = lapse.dydx[eh->z_ind][0][0]; //Xloc / (2.*r*r);//
	const double dXdr = X.dydx[eh->z_ind][0][0]; //-Xloc*Xloc*Xloc / (2.*r*r);//

	double tmp;
	double* xup = eh->p.xup;
	Tuple<double,40>& ch = eh->christoffel.data;

	// spatial parts
	unsigned offset;
	for(int a=0; a<3; a++){
		offset = 10*a;
		ch[offset+itt] = alpha * dadr / (r*Xloc*Xloc) * xup[a];
		ch[offset+ixt] = 0;
		ch[offset+iyt] = 0;
		ch[offset+izt] = 0;

		tmp = (1. - Xloc*Xloc + r*Xloc*dXdr) / (r*r*r*Xloc*Xloc) * xup[a]/r;
		ch[offset+ixx] = xup[0]*xup[0] * tmp;
		ch[offset+iyy] = xup[1]*xup[1] * tmp;
		ch[offset+izz] = xup[2]*xup[2] * tmp;
		ch[offset+ixy] = xup[0]*xup[1] * tmp;
		ch[offset+ixz] = xup[0]*xup[2] * tmp;
		ch[offset+iyz] = xup[1]*xup[2] * tmp;

		tmp = -(1.-Xloc*Xloc) / (r*Xloc*Xloc) * xup[a]/r;
		ch[offset+ixx] += tmp;
		ch[offset+iyy] += tmp;
		ch[offset+izz] += tmp;
	}
	// time part
	offset = 30;
	tmp = dadr / (r * alpha);
	ch[offset+itt] = 0;
	ch[offset+ixt] = xup[0] * tmp;
	ch[offset+iyt] = xup[1] * tmp;
	ch[offset+izt] = xup[2] * tmp;
	ch[offset+ixx] = 0;
	ch[offset+iyy] = 0;
	ch[offset+izz] = 0;
	ch[offset+ixy] = 0;
	ch[offset+ixz] = 0;
	ch[offset+iyz] = 0;


	for(unsigned i=0; i<40; i++) PRINT_ASSERT(ch[i],==,ch[i]);
}

void Grid1DSphere::axis_vector(vector<Axis>& axes) const{
	axes = vector<Axis>({rAxis});
}
double Grid1DSphere::zone_lorentz_factor(const int z_ind) const{
	double vdotv = vr[z_ind]*vr[z_ind] * X[z_ind] / (pc::c*pc::c);
	return 1. / sqrt(1.-vdotv);
}
void Grid1DSphere::interpolate_shift(const double xup[4], double betaup[3], const unsigned dir_ind[NDIMS]) const{ // default Minkowski
	betaup[0] = 0;
	betaup[1] = 0;
	betaup[2] = 0;
}
void Grid1DSphere::grid_coordinates(const double xup[3], double coords[NDIMS]) const{
	coords[0] = radius(xup);
}
