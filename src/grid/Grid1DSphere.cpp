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
	grid_type = "Grid1DSphere";
	reflect_outer = 0;
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
}

void Grid1DSphere::read_nagakura_model(Lua* lua){
	// verbocity
	int MPI_myID;
	MPI_Comm_rank( MPI_COMM_WORLD, &MPI_myID );
	const int rank0 = (MPI_myID == 0);
	vector<double> bintops;
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
	while(rgrid_file >> trash >> tmp){
		if(bintops.size()>0) PRINT_ASSERT(tmp,>,bintops[bintops.size()-1]);
		else PRINT_ASSERT(tmp,>,minval);
		bintops.push_back(tmp);
	}
	r_out.init(bintops,minval);

	// write grid properties
	if(rank0)
	  cout << "#   nr=" << r_out.size() << "\trmin=" << r_out.min << "\trmax=" << r_out[r_out.size()-1] << endl;

	// read the fluid properties
	z.resize(r_out.size());
	for(int z_ind=0; z_ind<r_out.size(); z_ind++){
		double trash;

		// read the contents of a single line
		infile >> trash; // r
		infile >> trash; // theta
		infile >> z[z_ind].rho; // g/ccm
		infile >> z[z_ind].Ye;
		infile >> z[z_ind].T; // MeV
		infile >> z[z_ind].u[0]; // cm/s
		infile >> z[z_ind].u[1]; // 1/s
		infile >> z[z_ind].u[2]; // 1/s

		// get rid of the rest of the line
		for(int k=9; k<=165; k++) infile >> trash;

		// convert units
		z[z_ind].u[1] *= r_out.center(z_ind);
		z[z_ind].u[2] *= r_out.center(z_ind);
		z[z_ind].T /= pc::k_MeV;

		// sanity checks
		PRINT_ASSERT(z[z_ind].rho,>=,0.0);
		PRINT_ASSERT(z[z_ind].T,>=,0.0);
		PRINT_ASSERT(z[z_ind].Ye,>=,0.0);
		PRINT_ASSERT(z[z_ind].Ye,<=,1.0);
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
	z.resize(n_zones);
	r_out.resize(n_zones);
	metric.resize(n_zones);

	// read zone properties
	vector<double> tmp_alpha = vector<double>(n_zones,0);
	vector<double> tmp_X = vector<double>(n_zones,0);
	infile >> r_out.min;
	PRINT_ASSERT(r_out.min,>=,0);
	for(int z_ind=0; z_ind<n_zones; z_ind++)
	{
		double alpha, x;
		infile >> r_out[z_ind];
		infile >> z[z_ind].rho;
		infile >> z[z_ind].T;
		infile >> z[z_ind].Ye;
		z[z_ind].H_vis = 0;
		infile >> z[z_ind].u[0];
		infile >> z[z_ind].u[1];
		infile >> z[z_ind].u[2];
		infile >> tmp_alpha[z_ind];
		infile >> tmp_X[z_ind];

		PRINT_ASSERT(r_out[z_ind],>,(z_ind==0 ? r_out.min : r_out[z_ind-1]));
		PRINT_ASSERT(z[z_ind].rho,>=,0);
		PRINT_ASSERT(z[z_ind].T,>=,0);
		PRINT_ASSERT(z[z_ind].Ye,>=,0);
		PRINT_ASSERT(z[z_ind].Ye,<=,1.0);
		PRINT_ASSERT(z[z_ind].u[0]*z[z_ind].u[0] + z[z_ind].u[1]*z[z_ind].u[1] + z[z_ind].u[2]*z[z_ind].u[2],<,pc::c*pc::c);
		PRINT_ASSERT(tmp_alpha[z_ind],<=,1.0);
		PRINT_ASSERT(tmp_X[z_ind],>=,1.0);
	}

	// set the christoffel symbol coefficients
	vector<double> tmp_dadr = vector<double>(n_zones,0);
	vector<double> tmp_dXdr = vector<double>(n_zones,0);
	if(do_GR) for(int z_ind=0; z_ind<metric.size(); z_ind++){
		double rL=0, r0=0, rR=0;
		double aL=0, a0=0, aR=0;
		double XL=0, X0=0, XR=0;
		double drL=0, drR=0;
		double daL=0, daR=0;
		double dXL=0, dXR=0;

		zone_coordinates(z_ind, &r0, 1);
		a0 = tmp_alpha[z_ind];
		X0 = tmp_X[z_ind];
		PRINT_ASSERT(X0,>,0);
		PRINT_ASSERT(a0,>,0);
		if(z_ind>0){
			zone_coordinates(z_ind-1, &rL, 1);
			aL = tmp_alpha[z_ind-1];
			XL = tmp_X[z_ind-1];
			drL = r0-rL;
			daL = a0-aL;
			dXL = X0-XL;
		}
		if(z_ind < z.size()-1){
			zone_coordinates(z_ind+1, &rR, 1);
			aR = tmp_alpha[z_ind+1];
			XR = tmp_X[z_ind+1];
			drR = rR-r0;
			daR = aR-a0;
			dXR = XR-X0;
		}

		double dX_dr=0, da_dr=0;

		if(z_ind == 0){
			da_dr = daR / drR;
			dX_dr = dXR / drR;
		}
		else if(z_ind == z.size()-1){
			da_dr = daL / drL;
			dX_dr = dXL / drL;
		}
		else{
			da_dr = (drR*daL/drL + drL*daR/drR) / (drL+drR);
			dX_dr = (drR*dXL/drL + drL*dXR/drR) / (drL+drR);
		}

		tmp_dadr[z_ind] = da_dr;
		tmp_dXdr[z_ind] = dX_dr;
	}
	metric.set_metric(tmp_alpha,tmp_X,tmp_dadr,tmp_dXdr);

	infile.close();
}


//------------------------------------------------------------
// Return the zone index containing the position x
//------------------------------------------------------------
int Grid1DSphere::zone_index(const double x[3]) const
{
	PRINT_ASSERT(z.size(),>,0);
	double r = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
	PRINT_ASSERT(r,>=,0);

	// check if off the boundaries
	if(r < r_out.min             ) return -1;
	if(r >= r_out[r_out.size()-1] ) return -2;

	// find in zone array using stl algorithm upper_bound and subtracting iterators
	int z_ind = r_out.locate(r);
	PRINT_ASSERT(z_ind,>=,0);
	if(z_ind>=(int)z.size()){
		cout << z_ind << endl;
		cout << z.size() << endl;
		cout << r << endl;
		cout << r_out[0] << endl;
	}
	PRINT_ASSERT(z_ind,<,(int)z.size());
	return z_ind;
}


//------------------------------------------------------------
// return volume of zone z_ind
//------------------------------------------------------------
double  Grid1DSphere::zone_lab_3volume(const int z_ind) const
{
	PRINT_ASSERT(z_ind,>=,0);
	PRINT_ASSERT(z_ind,<,(int)z.size());
	double r0 = (z_ind==0 ? r_out.min : r_out[z_ind-1]);
	double vol = 4.0*pc::pi/3.0*(r_out[z_ind]*r_out[z_ind]*r_out[z_ind] - r0*r0*r0);
	if(do_GR) vol *= metric.X[z_ind];
	PRINT_ASSERT(vol,>=,0);
	return vol;
}

//------------------------------------------------------------
// return length of zone
//------------------------------------------------------------
double  Grid1DSphere::zone_min_length(const int z_ind) const
{
	PRINT_ASSERT(z_ind,>=,0);
	PRINT_ASSERT(z_ind,<,(int)z.size());
	double r0 = (z_ind==0 ? r_out.min : r_out[z_ind-1]);
	double min_len = r_out[z_ind] - r0;
	PRINT_ASSERT(min_len,>=,0);
	return min_len;
}


// ------------------------------------------------------------
// find the coordinates of the zone in geometrical coordinates
// ------------------------------------------------------------
void Grid1DSphere::zone_coordinates(const int z_ind, double r[1], const int rsize) const{
	PRINT_ASSERT(z_ind,>=,0);
	PRINT_ASSERT(z_ind,<,(int)z.size());
	PRINT_ASSERT(rsize,==,(int)dimensionality());
	r[0] = 0.5*(r_out[z_ind]+r_out.bottom(z_ind));
	PRINT_ASSERT(r[0],>,0);
	PRINT_ASSERT(r[0],<,r_out[r_out.size()-1]);
}


//-------------------------------------------
// get directional indices from zone index
//-------------------------------------------
void Grid1DSphere::zone_directional_indices(const int z_ind, int dir_ind[1], const int size) const
{
	PRINT_ASSERT(z_ind,>=,0);
	PRINT_ASSERT(z_ind,<,(int)z.size());
	PRINT_ASSERT(size,==,(int)dimensionality());
	dir_ind[0] = z_ind;
}


//------------------------------------------------------------
// sample a random position within the spherical shell
//------------------------------------------------------------
void Grid1DSphere::sample_in_zone(const int z_ind, ThreadRNG* rangen, double x[3]) const
{
	PRINT_ASSERT(z_ind,>=,0);
	PRINT_ASSERT(z_ind,<,(int)z.size());

	double rand[3];
	rand[0] = rangen->uniform();
	rand[1] = rangen->uniform();
	rand[2] = rangen->uniform();

	// inner and outer radii of shell
	double r0 = (z_ind==0 ? r_out.min : r_out[z_ind-1]);
	double r1 = r_out[z_ind];

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
void Grid1DSphere::interpolate_fluid_velocity(const double x[3], double v[3], int z_ind) const
{
	if(z_ind < 0) z_ind = zone_index(x);
	PRINT_ASSERT(z_ind,>=,0);
	PRINT_ASSERT(z_ind,<,(int)z.size());

	// radius in zone
	double r = sqrt(dot_Minkowski<3>(x,x));

	// assuming radial velocity (may want to interpolate here)
	// (the other two components are ignored and mean nothing)
	v[0] = x[0]/r*z[z_ind].u[0];
	v[1] = x[1]/r*z[z_ind].u[0];
	v[2] = x[2]/r*z[z_ind].u[0];

	// check for pathological case
	if (r == 0)
	{
		v[0] = 0;
		v[1] = 0;
		v[2] = 0;
	}

	PRINT_ASSERT(dot_Minkowski<3>(v,v),<=,pc::c*pc::c);
}


//------------------------------------------------------------
// Write the grid information out to a file
//------------------------------------------------------------
void Grid1DSphere::write_rays(const int iw) const
{
	// this is a 1D grid, so the function is exactly the same
	// as write_zones
}


//------------------------------------------------------------
// Reflect off symmetry axis
//------------------------------------------------------------
void Grid1DSphere::symmetry_boundaries(LorentzHelper *lh) const{
	// reflect from outer boundary
	if(reflect_outer && radius(lh->p_xup())>r_out[r_out.size()-1]){
		const Particle *p = lh->particle_readonly(lab);

		double r0 = (r_out.size()>1 ? r_out[r_out.size()-2] : r_out.min);
		double rmax = r_out[r_out.size()-1];
		double dr = rmax - r0;
		double R = radius(p->xup);
		PRINT_ASSERT( fabs(R - r_out[r_out.size()-1]),<,TINY*dr);

		double kr;
		for(int i=0; i<3; i++) kr += p->xup[i]/R * p->kup[i];

		// invert the radial component of the velocity
		double knew[3];
		knew[3] = lh->p_kup(lab)[3];
		knew[0] -= 2.*kr * p->xup[0]/R;
		knew[1] -= 2.*kr * p->xup[1]/R;
		knew[2] -= 2.*kr * p->xup[2]/R;
		normalize_null(knew,p->xup);

		lh->set_p_kup<lab>(knew,4);

		// put the particle just inside the boundary
		double newR = rmax - TINY*dr;
		double x[4];
		x[0] = p->xup[0]/R*newR;
		x[1] = p->xup[1]/R*newR;
		x[2] = p->xup[2]/R*newR;
		x[3] = p->xup[3];
		lh->set_p_xup(x,4);

		// must be inside the boundary, or will get flagged as escaped
		PRINT_ASSERT(zone_index(x),>=,0);
	}
}

double Grid1DSphere::zone_cell_dist(const double x_up[3], const int z_ind) const{
	double r = sqrt(dot_Minkowski<3>(x_up,x_up));
	PRINT_ASSERT(r,<=,r_out[z_ind]);
	PRINT_ASSERT(r,>=,r_out.bottom(z_ind));

	double drL = r - r_out.bottom(z_ind);
	double drR = r_out[z_ind] - r;

	return min(drL, drR);
}

double Grid1DSphere::zone_radius(const int z_ind) const{
	PRINT_ASSERT(z_ind,>=,0);
	PRINT_ASSERT(z_ind,<,(int)z.size());
	return r_out[z_ind];
}

//-----------------------------
// Dimensions of the grid
//-----------------------------
void Grid1DSphere::dims(hsize_t dims[1], const int size) const{
	PRINT_ASSERT(size,==,(int)dimensionality());
	dims[0] = r_out.size();
}

//----------------------------------------------------
// Write the coordinates of the grid points to the hdf5 file
//----------------------------------------------------
void Grid1DSphere::write_hdf5_coordinates(H5::H5File file) const
{
	// useful quantities
	H5::DataSet dataset;
	H5::DataSpace dataspace;
	float tmp[r_out.size()+1];

	// get dimensions
	hsize_t coord_dims[1];
	dims(coord_dims,1);
	for(unsigned i=0; i<1; i++) coord_dims[i]++; //make room for min value

	// write coordinates
	dataspace = H5::DataSpace(1,&coord_dims[0]);
	dataset = file.createDataSet("r(cm)",H5::PredType::IEEE_F32LE,dataspace);
	tmp[0] = r_out.min;
	for(unsigned i=1; i<r_out.size()+1; i++) tmp[i] = r_out[i-1];
	dataset.write(&tmp[0],H5::PredType::IEEE_F32LE);
	dataset.close();
}

double Grid1DSphere::lapse(const double xup[4], int z_ind) const{
	double r = radius(xup);
	if(z_ind<0) z_ind = zone_index(xup);
	return metric.get_alpha(z_ind, r, r_out);
}

void Grid1DSphere::shiftup(double betaup[4], const double xup[4], int z_ind) const{
	for(int i=0; i<4; i++) betaup[i] = 0;
}

void Grid1DSphere::g3_down(const double xup[4], double gproj[4][4], int z_ind) const{
	const double r = radius(xup);
	if(z_ind<0) z_ind = zone_index(xup);
	const double X = metric.get_X(z_ind,r, r_out);
	for(int i=0; i<3; i++){
		for(int j=0; j<3; j++) gproj[i][j] = xup[i] * xup[j] * (X*X-1.0) / (r*r);
		gproj[i][i] += 1.0;
	}
}

void Grid1DSphere::connection_coefficients(const double xup[4], double gamma[4][4][4], int z_ind) const{
	if(z_ind<0) z_ind = zone_index(xup);
	const double r = radius(xup);
	const double X = metric.get_X(z_ind, r, r_out); // 1.0/alpha;
	const double alpha = metric.get_alpha(z_ind, r, r_out); // sqrt(1.0 - 1.0/r);
	const double dadr = metric.get_dadr(z_ind); // 0.5/(r*r) * X;
	const double dXdr = metric.get_dXdr(z_ind); // -0.5/(r*r) * X*X*X;

	// t mu nu
	gamma[3][3][3] = 0;
	for(int i=0; i<3; i++){
		for(int j=0; j<3; j++) gamma[3][i][j] = 0;
		gamma[3][i][3] = gamma[3][3][i] = dadr / (r*alpha) * xup[i];
	}

	for(int a=0; a<4; a++) for(int mu=0; mu<4; mu++) for(int nu=0; nu<4; nu++){
		double result=0.;
		if(a==3){
			if(mu==3 and nu==3) result=0;
			else if(nu==3) result =  dadr / (r*alpha) * xup[mu];
			else if(mu==3) result =  dadr / (r*alpha) * xup[nu];
			else result = 0;
		}

		else{ // a == 0-2
			if(mu==3 and nu==3) result = alpha * dadr / (r*X*X) * xup[a];
			else if(mu==3 or nu==3) result = 0;
			else{
				result = xup[mu]*xup[nu] / (r*r*r*X*X) * (1.0 - X*X + r*X*dXdr);
				if(mu==nu) result -= (1.-X*X) / (r*X*X);
				result *= xup[a] / r;
			}
		}

		gamma[a][mu][nu] = result;
	}
}

double Grid1DSphere::zone_lapse(const int z_ind) const{
	return metric.alpha[z_ind];
}
