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

	reflect_outer = lua->scalar<int>("reflect_outer");
}

void Grid1DSphere::write_child_zones(H5::H5File file){
	vr.write_HDF5(file,"vr(cm|s)");
	X.write_HDF5(file,"X");
}

void Grid1DSphere::read_nagakura_model(Lua* lua){
	// verbocity
	int MPI_myID;
	MPI_Comm_rank( MPI_COMM_WORLD, &MPI_myID );
	int rank0 = (MPI_myID == 0);
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
	xAxes[0] = Axis(minval, bintops, binmid);
	vector<Axis> axes = {xAxes[0]};
	vr.set_axes(axes);
	lapse.set_axes(axes);
	X.set_axes(axes);
	rho.set_axes(axes);
	T.set_axes(axes);
	Ye.set_axes(axes);
	H_vis.set_axes(axes);

	// write grid properties
	if(rank0)
	  cout << "#   nr=" << xAxes[0].size() << "\trmin=" << xAxes[0].min << "\trmax=" << xAxes[0].top[xAxes[0].size()-1] << endl;

	// read the fluid properties
	for(size_t z_ind=0; z_ind<xAxes[0].size(); z_ind++){
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
		T[z_ind] /= pc::k_MeV;

		// GR variables
		lapse[z_ind] = 1.0;
		X[z_ind] = 1.0;

		// sanity checks
		PRINT_ASSERT(rho[z_ind],>=,0.0);
		PRINT_ASSERT(T[z_ind],>=,0.0);
		PRINT_ASSERT(Ye[z_ind],>=,0.0);
		PRINT_ASSERT(Ye[z_ind],<=,1.0);
		PRINT_ASSERT(vr[z_ind],<,pc::c);
	}
}

void Grid1DSphere::read_custom_model(Lua* lua){
	// verbocity
	int MPI_myID;
	MPI_Comm_rank( MPI_COMM_WORLD, &MPI_myID );
	int rank0 = (MPI_myID == 0);
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
	xAxes[0] = Axis(rmin, rtop, rmid);
	vector<Axis> axes = {xAxes[0]};
	vr.set_axes(axes);
	lapse.set_axes(axes);
	X.set_axes(axes);
	rho.set_axes(axes);
	T.set_axes(axes);
	Ye.set_axes(axes);
	H_vis.set_axes(axes);

	for(size_t z_ind=0; z_ind<vr.size(); z_ind++){
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
int Grid1DSphere::zone_index(const Tuple<double,4>& x) const
{
	PRINT_ASSERT(rho.size(),>,0);
	double r = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
	PRINT_ASSERT(r,>=,0);

	// check if off the boundaries
	if(r < xAxes[0].min             ) return -1;
	if(r >= xAxes[0].top[xAxes[0].size()-1] ) return -1;

	// find in zone array using stl algorithm upper_bound and subtracting iterators
	int z_ind = xAxes[0].bin(r);
	PRINT_ASSERT(z_ind,>=,0);
	PRINT_ASSERT(z_ind,<,(int)rho.size());
	return z_ind;
}


//------------------------------------------------------------
// return volume of zone z_ind
//------------------------------------------------------------
double  Grid1DSphere::zone_lab_3volume(int z_ind) const
{
	PRINT_ASSERT(z_ind,>=,0);
	PRINT_ASSERT(z_ind,<,(int)rho.size());
	double r0 = (z_ind==0 ? xAxes[0].min : xAxes[0].top[z_ind-1]);
	double vol = 4.0*pc::pi/3.0*( pow(xAxes[0].top[z_ind],3) - pow(r0,3) );
	if(DO_GR) vol *= X[z_ind];
	PRINT_ASSERT(vol,>=,0);
	return vol;
}

//------------------------------------------------------------
// return length of zone
//------------------------------------------------------------
double  Grid1DSphere::zone_min_length(int z_ind) const
{
	PRINT_ASSERT(z_ind,>=,0);
	PRINT_ASSERT(z_ind,<,(int)rho.size());
	double r0 = (z_ind==0 ? xAxes[0].min : xAxes[0].top[z_ind-1]);
	double min_len = xAxes[0].top[z_ind] - r0;
	PRINT_ASSERT(min_len,>=,0);
	return min_len;
}


// ------------------------------------------------------------
// find the coordinates of the zone in geometrical coordinates
// ------------------------------------------------------------
Tuple<double,NDIMS> Grid1DSphere::zone_coordinates(int z_ind) const{
	PRINT_ASSERT(z_ind,>=,0);
	PRINT_ASSERT(z_ind,<,(int)rho.size());
	double r = 0.5*(xAxes[0].top[z_ind]+xAxes[0].bottom(z_ind));
	PRINT_ASSERT(r,>,0);
	PRINT_ASSERT(r,<,xAxes[0].top[xAxes[0].size()-1]);

	return Tuple<double,NDIMS>(r);
}


//-------------------------------------------
// get directional indices from zone index
//-------------------------------------------
Tuple<size_t,NDIMS> Grid1DSphere::zone_directional_indices(int z_ind) const
{
	PRINT_ASSERT(z_ind,>=,0);
	PRINT_ASSERT(z_ind,<,(int)rho.size());
	return Tuple<size_t,NDIMS>(z_ind);
}


//------------------------------------------------------------
// sample a random position within the spherical shell
//------------------------------------------------------------
Tuple<double,4> Grid1DSphere::sample_in_zone(int z_ind, ThreadRNG* rangen) const
{
	PRINT_ASSERT(z_ind,>=,0);
	PRINT_ASSERT(z_ind,<,(int)rho.size());

	double rand[3];
	rand[0] = rangen->uniform();
	rand[1] = rangen->uniform();
	rand[2] = rangen->uniform();

	// inner and outer radii of shell
	double r0 = (z_ind==0 ? xAxes[0].min : xAxes[0].top[z_ind-1]);
	double r1 = xAxes[0].top[z_ind];

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
	Tuple<double,4> x;
	x[0] = radius*sin_theta*cos(phi);
	x[1] = radius*sin_theta*sin(phi);
	x[2] = radius*mu;
	return x;
}


//------------------------------------------------------------
// get the velocity vector 
//------------------------------------------------------------
Tuple<double,3> Grid1DSphere::interpolate_fluid_velocity(const EinsteinHelper& eh) const
{
	// radius in zone
	double r = radius(eh.xup);

	if(r==0)
		return Tuple<double,3>(0);
	else{
		Tuple<double,3> x3;
		for(size_t i=0; i<3; i++) x3[i] = eh.xup[i];
		return x3 / r * vr.interpolate(eh.icube_vol);
	}
}



//------------------------------------------------------------
// Reflect off symmetry axis
//------------------------------------------------------------
bool reject_direction_inward(const EinsteinHelper* eh, ThreadRNG* rangen){
	double xdotx = eh->g.dot<3>(eh->xup, eh->xup);
	double kdotk = eh->g.dot<3>(eh->kup, eh->kup);
	double xdotk = eh->g.dot<3>(eh->xup, eh->kup);
	double costheta = xdotk / sqrt(xdotx * kdotk);
	return -costheta < rangen->uniform();
}
void Grid1DSphere::symmetry_boundaries(EinsteinHelper *eh) const{
	// reflect from outer boundary
	double R = radius(eh->xup);
	if(reflect_outer && R>xAxes[0].top[xAxes[0].size()-1]){
		double r0 = (xAxes[0].size()>1 ? xAxes[0].top[xAxes[0].size()-2] : xAxes[0].min);
		double rmax = xAxes[0].max();
		double dr = rmax - r0;

		double kr = 0;
		for(int i=0; i<3; i++) kr += eh->xup[i]/R * eh->kup[i];

		// give the particle an inward-moving direction
		Tuple<double,4> kup_tet;
		do{
			sim->isotropic_kup_tet(eh->nu(),kup_tet,&sim->rangen);
			eh->set_kup_tet(kup_tet);
		} while(reject_direction_inward(eh,&sim->rangen));
		//eh->g.normalize_null_preservedownt(eh->kup);

		// put the particle just inside the boundary
		double newR = rmax - TINY*dr;
		for(size_t i=0; i<3; i++)	eh->xup[i] *= newR/R;

		// must be inside the boundary, or will get flagged as escaped
		PRINT_ASSERT(zone_index(eh->xup),>=,0);
	}
}

double Grid1DSphere::d_boundary(const EinsteinHelper& eh) const{
	double r = radius(eh.xup);
	PRINT_ASSERT(r,<=,xAxes[0].top[eh.z_ind]);
	PRINT_ASSERT(r,>=,xAxes[0].bottom(eh.z_ind));

	// get component of k in the radial direction
	double kr = eh.g.dot<4>(eh.e[2],eh.kup);

	double dlambda = INFINITY;
	if(kr>0) dlambda = (xAxes[0].top[eh.z_ind]    - r) / kr;
	if(kr<0) dlambda = (xAxes[0].bottom(eh.z_ind) - r) / kr;

	double ds_com = dlambda * eh.kup_tet[3];
	PRINT_ASSERT(ds_com,>=,0);
	return ds_com;
}
double Grid1DSphere::d_randomwalk(const EinsteinHelper& eh) const{
	double R=INFINITY;
	double D = pc::c / (3.*eh.scatopac);

	Tuple<double,4> ktest;
	for(size_t i=0; i<3; i++) ktest[i] = eh.xup[i];
	ktest[3] = 0;
	const double r = radius(eh.xup);
	const double kr = r;
	const double ur = radius(eh.u);

	for(int sgn=1; sgn>0; sgn*=-1){
		// get a null test vector
		for(size_t i=0; i<3; i++) ktest[i] *= sgn;
		eh.g.normalize_null_changeupt(ktest);

		// get the time component of the tetrad test vector
		double kup_tet_t = -eh.g.dot<4>(ktest,eh.u);

		// get the min distance from the boundary in direction i. Negative if moving left
		double drlab=0;
		if(sgn>0) drlab = xAxes[0].top[eh.dir_ind[0]] - r;
		if(sgn<0) drlab = xAxes[0].bottom(eh.dir_ind[0]) - r;

		R = min(R, sim->R_randomwalk(kr/kup_tet_t, ktest[3]/kup_tet_t, ur, drlab, D));
	}

	PRINT_ASSERT(R,>=,0);
	PRINT_ASSERT(R,<,INFINITY);
	return R;
}

double Grid1DSphere::zone_radius(int z_ind) const{
	PRINT_ASSERT(z_ind,>=,0);
	PRINT_ASSERT(z_ind,<,(int)rho.size());
	return xAxes[0].top[z_ind];
}

//-----------------------------
// Dimensions of the grid
//-----------------------------
Tuple<hsize_t,NDIMS> Grid1DSphere::dims() const{
	return Tuple<hsize_t,NDIMS>(xAxes[0].size());
}


Tuple<double,6> Grid1DSphere::interpolate_3metric(const EinsteinHelper& eh) const{
	const double r = radius(eh.xup);
	const double Xloc = X.interpolate(eh.icube_vol);//1./sqrt(1.-1./r);//
	double tmp = (Xloc*Xloc-1.0) / (r*r);

	Tuple<double,6> data;
	data[ixx] = eh.xup[0]*eh.xup[0] * tmp;
	data[iyy] = eh.xup[1]*eh.xup[1] * tmp;
	data[izz] = eh.xup[2]*eh.xup[2] * tmp;
	data[ixy] = eh.xup[0]*eh.xup[1] * tmp;
	data[ixz] = eh.xup[0]*eh.xup[2] * tmp;
	data[iyz] = eh.xup[1]*eh.xup[2] * tmp;

	data[ixx] += 1.0;
	data[iyy] += 1.0;
	data[izz] += 1.0;

	return data;
}

Tuple<double,4> Grid1DSphere::dk_dlambda(const EinsteinHelper& eh) const{
	const double r = radius(eh.xup);
	const double alpha = lapse.interpolate(eh.icube_vol); //sqrt(1.-1./r); //
	const double Xloc  = X.interpolate(eh.icube_vol); //1./alpha; //
	const double dadr  = lapse.interpolate_slopes(eh.icube_vol)[0]; //Xloc / (2.*r*r);//
	const double dXdr  = X.interpolate_slopes(eh.icube_vol)[0]; //-Xloc*Xloc*Xloc / (2.*r*r);//

	double tmp;
	Christoffel ch;
	ch.data = 0;

	// spatial parts
	for(int a=0; a<3; a++){
		ch.data[Christoffel::index(a,3,3)] = alpha * dadr / (r*Xloc*Xloc) * eh.xup[a];

		tmp = (1. - Xloc*Xloc + r*Xloc*dXdr) / (r*r*r*Xloc*Xloc) * eh.xup[a]/r;
		for(size_t i=0; i<3; i++) for(size_t j=i; j<3; j++)
			ch.data[Christoffel::index(a,i,j)] = eh.xup[i]*eh.xup[j] * tmp;

		tmp = -(1.-Xloc*Xloc) / (r*Xloc*Xloc) * eh.xup[a]/r;
		for(size_t i=0; i<3; i++)
			ch.data[Christoffel::index(a,i,i)] += tmp;
	}
	// time part
	tmp = dadr / (r * alpha);
	for(size_t i=0; i<3; i++)
		ch.data[Christoffel::index(3,3,i)] = eh.xup[i] * tmp;

	for(size_t i=0; i<40; i++) PRINT_ASSERT(ch.data[i],==,ch.data[i]);

	return ch.contract2(eh.kup)*(-1);
}

double Grid1DSphere::zone_lorentz_factor(int z_ind) const{
	double vdotv = vr[z_ind]*vr[z_ind] * X[z_ind] / (pc::c*pc::c);
	return 1. / sqrt(1.-vdotv);
}
Tuple<double,3> Grid1DSphere::interpolate_shift(const EinsteinHelper&) const{ // default Minkowski
	return Tuple<double,3>(0);
}
void Grid1DSphere::grid_coordinates(const Tuple<double,4>& xup, double coords[NDIMS]) const{
	coords[0] = radius(xup);
}
