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

#ifndef _POLAR_SPECTRUM_ARRAY_H
#define _POLAR_SPECTRUM_ARRAY_H 1

#include "SpectrumArray.h"
#include "H5Cpp.h"
#include "Axis.h"
#include <mpi.h>
#include <sstream>
#include <fstream>
#include "global_options.h"
#include "MultiDArray.h"

using namespace std;
namespace pc = physical_constants;

template<unsigned ndims_spatial>
class PolarSpectrumArray : public SpectrumArray {

private:

	// bin arrays
	// values represent bin upper walls (the single locate_array.min value is the leftmost wall)
	// underflow is combined into leftmost bin (right of the locate_array.min)
	// overflow is combined into the rightmost bin (left of locate_array[size-1])


public:

	ScalarMultiDArray<ndims_spatial+3> data;
	unsigned phiGridIndex, nuGridIndex, muGridIndex;
	unsigned nphi, nnu, nmu;



	unsigned size() const{
		return data.size();
	}

	//--------------------------------------------------------------
	// Initialization and Allocation
	//--------------------------------------------------------------
	void init(const vector<Axis>& spatial_axes,
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
		data.set_axes(axes);
		data.wipe();
	}


	//--------------------------------------------------------------
	// Initialization and Allocation
	//--------------------------------------------------------------
	void init(const vector<Axis>& spatial_axes, const Axis& wg,	const Axis& mg, const Axis& pg){
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
		data.set_axes(axes);
		data.wipe();
	}


	//--------------------------------------------------------------
	// Functional procedure: Wipe
	//--------------------------------------------------------------
	void wipe(){
		data.wipe();
	}


	//--------------------------------------------------------------
	// count a particle
	////--------------------------------------------------------------
	void count(const EinsteinHelper* eh, const double E){
		PRINT_ASSERT(E,>=,0);
		PRINT_ASSERT(eh->kup_tet[3],>=,0);
		const double tiny = 1e-8;

		unsigned indices[data.Ndims()];
		for(int i=0; i<ndims_spatial; i++) indices[i] = eh->dir_ind[i];
		indices[nuGridIndex] = eh->dir_ind[NDIMS];
		
		double mu = eh->kup_tet[2] / eh->kup_tet[3];
		mu = max(-1.0+tiny,mu);
		mu = min( 1.0-tiny,mu);
		int mu_bin = data.axes[muGridIndex].bin(mu);
		mu_bin = max(mu_bin, 0);
		mu_bin = min(mu_bin, (int)data.axes[muGridIndex].size()-1);
		indices[muGridIndex] = mu_bin;

		double phi = atan2(eh->kup_tet[1],eh->kup_tet[0]);  // projection into x-y plane
		if(phi< -pc::pi) phi += 2.0*pc::pi;
		if(phi>= pc::pi) phi -= 2.0*pc::pi;
		int phi_bin = data.axes[phiGridIndex].bin(phi);
		phi_bin = max(phi_bin, 0);
		phi_bin = min((unsigned)phi_bin, data.axes[phiGridIndex].size()-1);
		indices[phiGridIndex] = phi_bin;

		data.add(indices, E);
	}

	void rescale(double r){
		for(unsigned i=0;i<data.size();i++) data.y0[i] *= r;
	}


	//--------------------------------------------------------------
	// MPI average the spectrum contents
	//--------------------------------------------------------------
	// only process 0 gets the reduced spectrum to print
	void MPI_average(){
		data.MPI_combine();
	}

	//--------------------------------------------------------------
	// Write data to specified location in an HDF5 file
	//--------------------------------------------------------------
	void write_hdf5_data(H5::H5File file, const string name) const{
		data.write_HDF5(file, name);
	}

	//--------------------------------------------------------------
	// Write distribution function coordinates to an HDF5 file
	//--------------------------------------------------------------
	void write_hdf5_coordinates(H5::H5File file, const string name) const{
		// useful quantities
		hsize_t dims[1];
		H5::DataSpace dataspace;
		H5::DataSet dataset;
		vector<float> tmp;

		data.axes[muGridIndex].write_HDF5(name+"_costheta_grid(lab)",file);
		data.axes[phiGridIndex].write_HDF5(name+"_phi_grid(radians,lab)",file);
	}
	
	void add_isotropic(const unsigned dir_ind[NDIMS+1], const double E){
		unsigned indices[data.Ndims()];
		for(int i=0; i<ndims_spatial; i++) indices[i] = dir_ind[i];
		indices[nuGridIndex] = dir_ind[ndims_spatial];
		indices[muGridIndex] = 0;
		indices[phiGridIndex] = 0;
		unsigned start = data.direct_index(indices);
		unsigned stop = start + nphi*nmu;
		double tmp = E / (double)(nphi*nmu);

		#pragma omp critical
		for(unsigned i=start; i<stop; i++) data[i] += tmp;
	}
};

#endif
