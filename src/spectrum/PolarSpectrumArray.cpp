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
#include <sstream>
#include <fstream>
#include "global_options.h"
#include "PolarSpectrumArray.h"
#include "Transport.h"

using namespace std;
namespace pc = physical_constants;

//--------------------------------------------------------------
// Initialization and Allocation
//--------------------------------------------------------------
void PolarSpectrumArray::init(const vector<Axis>& spatial_axes,
		const std::vector<double> w,
		const int n_mu, const int n_phi)
{
	vector<Axis> axes;

	// spatial axes
	for(int i=0; i<spatial_axes.size(); i++) axes.push_back(spatial_axes[i]);

	// frequency axis
	vector<double> top, mid;
	double min;
	double start = w[0];
	double stop  = w[1];
	double del   = w[2];
	unsigned ng = (stop-start)/del;
	axes.push_back(Axis(start, stop, ng));
	nuGridIndex = axes.size()-1;
	nnu = axes[axes.size()-1].size();

	// polar axes
	axes.push_back(Axis(-1,1,n_mu));
	muGridIndex = axes.size()-1;
	nmu = axes[axes.size()-1].size();
	axes.push_back(Axis(-pc::pi, pc::pi, n_phi));
	phiGridIndex = axes.size()-1;
	nphi = axes[axes.size()-1].size();

	// set up the data structure
	switch(spatial_axes.size()){
	case 0: data = new MultiDArray<0+3>(axes); break;
	case 1: data = new MultiDArray<1+3>(axes); break;
	case 2: data = new MultiDArray<2+3>(axes); break;
	case 3: data = new MultiDArray<3+3>(axes); break;
	}

	data->wipe();
}


//--------------------------------------------------------------
// Initialization and Allocation
//--------------------------------------------------------------
void PolarSpectrumArray::init(const vector<Axis>& spatial_axes, const Axis& wg,
		const Axis& mg, const Axis& pg)
{
	vector<Axis> axes;

	// spatial axes
	for(int i=0; i<spatial_axes.size(); i++) axes.push_back(spatial_axes[i]);

	axes.push_back(wg);
	nuGridIndex = axes.size()-1;
	nnu = wg.size();
	axes.push_back(mg);
	muGridIndex = axes.size()-1;
	nmu = mg.size();
	axes.push_back(pg);
	phiGridIndex = axes.size()-1;
	nphi = pg.size();

	// set up the data structure
	switch(spatial_axes.size()){
	case 0: data = new MultiDArray<0+3>(axes); break;
	case 1: data = new MultiDArray<1+3>(axes); break;
	case 2: data = new MultiDArray<2+3>(axes); break;
	case 3: data = new MultiDArray<3+3>(axes); break;
	}

	data->wipe();
}


//--------------------------------------------------------------
// Functional procedure: Wipe
//--------------------------------------------------------------
void PolarSpectrumArray::wipe()
{
	data->wipe();
}


//--------------------------------------------------------------
// count a particle
////--------------------------------------------------------------
void PolarSpectrumArray::count(const double D[3], const vector<unsigned>& dir_ind, const double nu, const double E)
{
	PRINT_ASSERT(E,>=,0);
	PRINT_ASSERT(nu,>=,0);
	PRINT_ASSERT(dir_ind.size(),==,data->Ndims()-3);
	const double tiny = 1e-8;

	unsigned indices[data->Ndims()];
	for(int i=0; i<dir_ind.size(); i++) indices[i] = dir_ind[i];

	double mu = D[2];
	mu = max(-1.0+tiny,mu);
	mu = min( 1.0-tiny,mu);
	int mu_bin = data->axes[muGridIndex].bin(mu);
	mu_bin = max(mu_bin, 0);
	mu_bin = min(mu_bin, data->axes[muGridIndex].size()-1);
	indices[muGridIndex] = mu_bin;

	double phi = atan2(D[1],D[0]);  // projection into x-y plane
	if(phi< -pc::pi) phi += 2.0*pc::pi;
	if(phi>= pc::pi) phi -= 2.0*pc::pi;
	int phi_bin = data->axes[phiGridIndex].bin(phi);
	phi_bin = max(phi_bin, 0);
	phi_bin = min(phi_bin, data->axes[phiGridIndex].size()-1);
	indices[phiGridIndex] = phi_bin;

	int nu_bin = data->axes[nuGridIndex].bin(nu);
	nu_bin = max(nu_bin, 0);
	nu_bin = min(nu_bin, data->axes[nuGridIndex].size()-1);
	indices[nuGridIndex] = nu_bin;

	data->add(indices, E);
}


void  PolarSpectrumArray::rescale(double r)
{
	for(unsigned i=0;i<data->size();i++) data->y0[i] *= r;
}


//--------------------------------------------------------------
// MPI average the spectrum contents
//--------------------------------------------------------------
// only process 0 gets the reduced spectrum to print
void PolarSpectrumArray::MPI_average()
{
	data->MPI_combine();
}

//--------------------------------------------------------------
// Write data to specified location in an HDF5 file
//--------------------------------------------------------------
void PolarSpectrumArray::write_hdf5_data(H5::H5File file, const string name) const
{
	data->write_HDF5(file, name);
}

//--------------------------------------------------------------
// Write distribution function coordinates to an HDF5 file
//--------------------------------------------------------------
void PolarSpectrumArray::write_hdf5_coordinates(H5::H5File file, const string name) const
{
	// useful quantities
	hsize_t dims[1];
	H5::DataSpace dataspace;
	H5::DataSet dataset;
	vector<float> tmp;

	data->axes[muGridIndex].write_HDF5(name+"_costheta_grid(lab)",file);
	data->axes[phiGridIndex].write_HDF5(name+"_phi_grid(radians,lab)",file);
}
