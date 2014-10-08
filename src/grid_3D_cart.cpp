#include "grid_3D_cart.h"
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iomanip>
#include "global_options.h"
#include "mpi.h"
#include "Lua.h"

//------------
// constructor
//------------
grid_3D_cart::grid_3D_cart(){
	nx=-MAX;ny=-MAX;nz=-MAX;
	dx=NaN;dy=NaN;dz=NaN;
	x0=NaN;y0=NaN;z0=NaN;
	vol=NaN;
	min_ds=NaN;
	reflect_x=-MAX;reflect_y=-MAX;reflect_z=-MAX;
}

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
	z.resize(n_zones,zone(dimensionality));
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
	z.resize(n_zones,zone(dimensionality));
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
	int i = floor((x[0]-x0)/dx);
	int j = floor((x[1]-y0)/dy);
	int k = floor((x[2]-z0)/dz);

	// check for off grid
	if ((i < 0)||(i > nx-1)) return -2;
	if ((j < 0)||(j > ny-1)) return -2;
	if ((k < 0)||(k > nz-1)) return -2;

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

	dir_ind.resize(dimensionality);
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
	return vol;
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
	vector<int> dir_ind;
	zone_directional_indices(z_ind,dir_ind);
	assert(dir_ind.size()==dimensionality);
	x[0] = x0 + ((double)dir_ind[0] + rand[0])*dx;
	x[1] = y0 + ((double)dir_ind[1] + rand[1])*dy;
	x[2] = z0 + ((double)dir_ind[2] + rand[2])*dz;
	assert(x[0] >= x0*(1.-tiny));
	assert(x[1] >= y0*(1.-tiny));
	assert(x[2] >= z0*(1.-tiny));
	assert(x[0] >= (x0+(double)dir_ind[0]*dx)*(1.+tiny));
	assert(x[1] >= (y0+(double)dir_ind[1]*dy)*(1.+tiny));
	assert(x[2] >= (z0+(double)dir_ind[2]*dz)*(1.+tiny));
	x[0] = max(x[0], x0);
	x[1] = max(x[1], y0);
	x[2] = max(x[2], z0);
	x[0] = min(x[0], x0+(double)dir_ind[0]*dx);
	x[1] = min(x[1], y0+(double)dir_ind[1]*dy);
	x[2] = min(x[2], z0+(double)dir_ind[2]*dz);
}


//------------------------------------------------------------
// return length of zone
//------------------------------------------------------------
double  grid_3D_cart::zone_min_length(const int z_ind) const
{
	assert(z_ind >= 0);
	assert(z_ind < (int)z.size());
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

	// may want to interpolate here
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
	r.resize(dimensionality);
	vector<int> dir_ind;
	zone_directional_indices(z_ind,dir_ind);
	r[0] = x0 + ((double)dir_ind[0]+0.5)*dx;
	r[1] = y0 + ((double)dir_ind[1]+0.5)*dy;
	r[2] = z0 + ((double)dir_ind[2]+0.5)*dz;
}


//------------------------------------------------------------
// Write the grid information out to a file
//------------------------------------------------------------
void grid_3D_cart::write_rays(const int iw) const
{
	int i,j,k;
	int z_ind;
	vector<double> r;

	ofstream outf;

	// XY-slice
	transport::open_file("slice_xy",iw,outf);
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
	transport::open_file("slice_xz",iw,outf);
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
	transport::open_file("slice_yz",iw,outf);
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
	transport::open_file("ray_x",iw,outf);
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
	transport::open_file("ray_y",iw,outf);
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
	transport::open_file("ray_z",iw,outf);
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
double grid_3D_cart::lab_dist_to_boundary(const particle *p) const{
	double xmax = x0 + nx*dx;
	double ymax = y0 + ny*dy;
	double zmax = z0 + nz*dz;

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
