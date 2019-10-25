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

template<size_t ndims_spatial>
class PolarSpectrumArray : public SpectrumArray {

private:

	// bin arrays
	// values represent bin upper walls (the single locate_array.min value is the leftmost wall)
	// underflow is combined into leftmost bin (right of the locate_array.min)
	// overflow is combined into the rightmost bin (left of locate_array[size-1])


public:

	ScalarMultiDArray<ATOMIC<double>,ndims_spatial+3> data;
	size_t phiGridIndex, nuGridIndex, muGridIndex;
	size_t nphi, nnu, nmu;

	size_t direct_index(const size_t dir_ind[ndims_spatial+3]) const{
		return data.direct_index(dir_ind);
	}
	double get(const size_t index) const{
		return data[index];
	}

	size_t size() const{
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
		size_t ng = (stop-start)/del;
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
		for(size_t i=0; i<spatial_axes.size(); i++) axes.push_back(spatial_axes[i]);

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
	void count_single(const Tuple<double,4>& kup_tet, const size_t dir_ind[NDIMS+1], const double E){
		PRINT_ASSERT(E,>=,0);
		PRINT_ASSERT(kup_tet[3],>=,0);

		size_t indices[data.Ndims()];
		for(size_t i=0; i<ndims_spatial; i++) indices[i] = dir_ind[i];
		indices[nuGridIndex] = dir_ind[NDIMS];
		
		double mu = kup_tet[2] / kup_tet[3];
		mu = max(-1.0+TINY,mu);
		mu = min( 1.0-TINY,mu);
		int mu_bin = data.axes[muGridIndex].bin(mu);
		mu_bin = max(mu_bin, 0);
		mu_bin = min(mu_bin, (int)data.axes[muGridIndex].size()-1);
		indices[muGridIndex] = mu_bin;

		double phi = atan2(kup_tet[1],kup_tet[0]);  // projection into x-y plane
		if(phi< -pc::pi) phi += 2.0*pc::pi;
		if(phi>= pc::pi) phi -= 2.0*pc::pi;
		int phi_bin = data.axes[phiGridIndex].bin(phi);
		phi_bin = max(phi_bin, 0);
		phi_bin = min(phi_bin, (int)data.axes[phiGridIndex].size()-1);
		indices[phiGridIndex] = phi_bin;

		data.add(indices, E);
	}

	void rescale(double r){
		for(size_t i=0;i<data.size();i++) data.y0[i] *= r;
	}
	void rescale_spatial_point(const size_t dir_ind[ndims_spatial], const double r){
		size_t all_indices[ndims_spatial+3];
		for(size_t i=0; i<ndims_spatial; i++) all_indices[i] = dir_ind[i];
		all_indices[ndims_spatial  ] = 0;
		all_indices[ndims_spatial+1] = 0;
		all_indices[ndims_spatial+2] = 0;
		size_t base_ind = data.direct_index(all_indices);
		size_t nbins = data.axes[ndims_spatial].size() * data.axes[ndims_spatial+1].size() * data.axes[ndims_spatial+2].size();
		for(size_t i=0; i<nbins; i++){
			data.y0[base_ind+i] *= r;
		}
	}

	//--------------------------------------------------------------
	// MPI scatter the spectrum contents.
	// Must rescale zone stop list to account for number of groups
	//--------------------------------------------------------------
	void mpi_sum_scatter(vector<size_t>& zone_stop_list){
		size_t nperzone = data.axes[nuGridIndex].size() * data.axes[phiGridIndex].size() * data.axes[muGridIndex].size();
		vector<size_t> stop_list = zone_stop_list;
		for(size_t i=0; i<stop_list.size(); i++) stop_list[i] *= nperzone;
		data.mpi_sum_scatter(stop_list);
	}
	void mpi_sum(){
		data.mpi_sum();
	}

	//--------------------------------------------------------------
	// Write data to specified location in an HDF5 file
	//--------------------------------------------------------------
	void write_hdf5_data(H5::H5File file, const string name) {
		data.write_HDF5(file, name);
	}
	void read_hdf5_data(H5::H5File file, const string name, const string axis_base) {
		vector<Axis> axes(ndims_spatial+3);
		for(hsize_t dir=0; dir<ndims_spatial; dir++)
			axes[dir].read_HDF5(string("/axes/x")+to_string(dir)+string("(cm)"),file);
		axes[nuGridIndex].read_HDF5("/axes/frequency(Hz)",file);
		axes[muGridIndex].read_HDF5(axis_base+"_costheta_grid(lab)",file);
		axes[phiGridIndex].read_HDF5(axis_base+"_phi_grid(radians,lab)",file);

		data.read_HDF5(file, name, axes);
	}

	//--------------------------------------------------------------
	// Write distribution function coordinates to an HDF5 file
	//--------------------------------------------------------------
	void write_hdf5_coordinates(H5::H5File file, const string name) const{
		// useful quantities
		H5::DataSpace dataspace;
		H5::DataSet dataset;
		vector<float> tmp;

		data.axes[muGridIndex].write_HDF5(name+"_costheta_grid(lab)",file);
		data.axes[phiGridIndex].write_HDF5(name+"_phi_grid(radians,lab)",file);
	}
	
	void add_isotropic_single(const size_t dir_ind[NDIMS+1], const double E){
		size_t indices[data.Ndims()];
		for(size_t i=0; i<ndims_spatial; i++) indices[i] = dir_ind[i];
		indices[nuGridIndex] = dir_ind[ndims_spatial];
		indices[muGridIndex] = 0;
		indices[phiGridIndex] = 0;
		size_t start = data.direct_index(indices);
		size_t stop = start + nphi*nmu;
		double tmp = E / (double)(nphi*nmu);

		#pragma omp critical
		for(size_t i=start; i<stop; i++) data[i] += tmp;
	}
	double total() const{
		double result=0;
		for(size_t i=0; i<data.size(); i++)
			result += data[i];
		return result;
	}

	static double cos_angle_between(const double mu1, const double mu2, const double phi1, const double phi2){
		PRINT_ASSERT(mu1,<=,1.0);
		PRINT_ASSERT(mu1,>=,-1.0);
		PRINT_ASSERT(mu2,<=,1.0);
		PRINT_ASSERT(mu2,>=,-1.0);
		double result = mu1*mu2 + sqrt((1.0-mu1*mu1)*(1.0-mu2*mu2))*cos(phi2-phi1); // Bruenn 1985 eq. A7

		// make sure numerical error doesn't push it beyond {-1,1}
		result = min(result, 1.0);
		result = max(result,-1.0);
		return result;
	}
	void annihilation_rate(
			const size_t dir_ind[NDIMS],       // directional indices for the zone we're getting the rate at
			const SpectrumArray* in_dist,  // erg/ccm (integrated over angular bin and energy bin)
			const vector< vector<vector<double> > >& phi, // cm^3/s [order][igin][igout]
			const size_t weight,
			Tuple<double,4>& fourforce) const{

		const PolarSpectrumArray<NDIMS>* nubar_dist = (PolarSpectrumArray<NDIMS>*)in_dist;
		PRINT_ASSERT(phi.size(),==,3);

		PRINT_ASSERT(size(),==,nubar_dist->size());
		PRINT_ASSERT(nnu,==,nubar_dist->nnu);
		PRINT_ASSERT(nmu,==,nubar_dist->nmu);
		PRINT_ASSERT(nphi,==,nubar_dist->nphi);

		const Axis* nuAxis = &(data.axes[nuGridIndex]);
		const Axis* muAxis = &(data.axes[muGridIndex]);
		const Axis* phiAxis = &(data.axes[phiGridIndex]);

		// calculate angle between distribution function angles beforehand
		double costheta[nmu][nphi][nmu][nphi];
		for(size_t mu=0; mu<nmu; mu++){
			for(size_t phi=0; phi<nphi; phi++){

				double avg_mu  = muAxis->mid[mu];
				double avg_phi = phiAxis->mid[phi];

				for(size_t mubar=0; mubar<nmu; mubar++){
					for(size_t phibar=0; phibar<nphi; phibar++){

						//size_t indexbar = nubar_dist->index(0,mubar,phibar);
						double avg_mubar  = muAxis->mid[mubar];
						double avg_phibar = phiAxis->mid[phibar];

						costheta[mu][phi][mubar][phibar] = cos_angle_between(avg_mu,avg_mubar,avg_phi,avg_phibar);
					}
				}
			}
		}

		// set up index arrays
		size_t index[ndims_spatial + 3];
		size_t indexbar[ndims_spatial + 3];
		for(size_t i=0; i<ndims_spatial; i++){
			index[i] = dir_ind[i];
			indexbar[i] = dir_ind[i];
		}

		// integrate over all bins
		for(size_t inu=0; inu<nnu; inu++){
			double avg_e = nuAxis->mid[inu]*pc::h; // erg
			index[nuGridIndex] = inu;
			for(size_t inubar=0; inubar<nnu; inubar++){
				double avg_ebar = nuAxis->mid[inubar]*pc::h; // erg
				indexbar[nuGridIndex] = inubar;

				double eebar = avg_e * avg_ebar;
				double phi2 = -1./5. * (phi[0][inu][inubar] + 3.*phi[1][inu][inubar]);

				// neutrino direction loops
				for(size_t imu=0; imu<nmu; imu++){
					index[muGridIndex] = imu;
					const double z = muAxis->mid[imu];
					const double sintheta = sqrt(1. - muAxis->mid[imu]*muAxis->mid[imu]);
					for(size_t iphi=0; iphi<nphi; iphi++){
						index[phiGridIndex] = iphi;
						size_t ind = direct_index(index);
						const double x = sintheta * cos(phiAxis->mid[iphi]);
						const double y = sintheta * sin(phiAxis->mid[iphi]);

						// antineutrino direction loops
						for(size_t imubar=0; imubar<nmu; imubar++){
							indexbar[nubar_dist->muGridIndex] = imubar;
							const double zbar = muAxis->mid[imubar];
							const double sinthetabar = sqrt(1. - muAxis->mid[imubar]*muAxis->mid[imubar]);
							for(size_t iphibar=0; iphibar<nphi; iphibar++){
								indexbar[nubar_dist->phiGridIndex] = iphibar;
								size_t indbar = nubar_dist->direct_index(indexbar);
								const double xbar = sinthetabar * cos(phiAxis->mid[iphibar]);
								const double ybar = sinthetabar * sin(phiAxis->mid[iphibar]);

								double cost = costheta[imu][iphi][imubar][iphibar];
								double nudist_edens    = get(ind); // erg/ccm
								double nubardist_edens = nubar_dist->get(indbar); // erg/ccm

								double Q = 0.5*phi[0][inu][inubar] +
										1.5*phi[1][inu][inubar] * cost +
										2.5*phi2 * 0.5*(3.*cost*cost-1.); // ccm/s
								Q *= nudist_edens * nubardist_edens / eebar; // #/ccm/s
								//PRINT_ASSERT(Q,>=,0);

								// angular distribution
								fourforce[0] += Q * (avg_e*x + avg_ebar*xbar); // erg/ccm/s
								fourforce[1] += Q * (avg_e*y + avg_ebar*ybar);
								fourforce[2] += Q * (avg_e*z + avg_ebar*zbar);
								fourforce[3] += Q * (avg_e   + avg_ebar     );

							} // phibar
						} // mubar
					} // phi
				} // mu
			} // nubar
		} // nu
		for(size_t i=0; i<4; i++) fourforce[i] /= weight;

		// sanity checks
		PRINT_ASSERT(fourforce[3],>=,0);
		for(size_t i=0; i<3; i++) PRINT_ASSERT(abs(fourforce[i]),<=,fourforce[3]);
	}

};

#endif
