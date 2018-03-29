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

#ifndef _MOMENT_SPECTRUM_ARRAY_H
#define _MOMENT_SPECTRUM_ARRAY_H 1

#include "SpectrumArray.h"
#include "H5Cpp.h"
#include <fstream>
#include <sstream>

using namespace std;
namespace pc = physical_constants;

template<unsigned ndims_spatial>
class MomentSpectrumArray : public SpectrumArray {

private:

	// helper variables
	static const unsigned index_range = 3; // index can be 0,1,2
	static const unsigned nranks = 4;
	static const unsigned nuGridIndex = ndims_spatial;
	const unsigned n_independent_elements[4] = {1,3,6,10}; // for a symmetric tensor of rank 0,1,2,3
	static const unsigned n_total_elements = 1+3+6+10;

	// bin arrays
	// values represent bin upper walls (the single locate_array.min value is the leftmost wall)
	// underflow is combined into leftmost bin (right of the locate_array.min)
	// overflow is combined into the rightmost bin (left of locate_array[size-1])
	MultiDArray<double,n_total_elements, ndims_spatial+1> data;

public:

	unsigned direct_index(const unsigned dir_ind[ndims_spatial+1]) const{
		return data.direct_index(dir_ind);
	}
	double getE(const unsigned ind) const{
		return data[ind][0];
	}
	double getF(const unsigned ind, const unsigned i) const{
		return data[ind][i+1];
	}
	double getP(const unsigned ind, const unsigned i, const unsigned j) const{
		switch((i+1)*(j+1)){
		case 1:
			return data[ind][4];
			break;
		case 2:
			return data[ind][5];
			break;
		case 3:
			return data[ind][6];
			break;
		case 4:
			return data[ind][7];
			break;
		case 6:
			return data[ind][8];
			break;
		case 9:
			return data[ind][9];
			break;
		}

		// should not get past switch statement
		assert(0);
		return NaN;
	}
	double getL(const unsigned ind, const unsigned i, const unsigned j, const unsigned k) const{
		switch((i+1)*(j+1)*(k+1)){
		case 1:
			return data[ind][10]; //Lxxx=data[index][10]
			break;
		case 2:
			return data[ind][11]; //Lxxy=data[index][11]
			break;
		case 3:
			return data[ind][12]; //Lxxz=data[index][12]
			break;
		case 4:
			return data[ind][13]; //Lxyy=data[index][13]
			break;
		case 6:
			return data[ind][14]; //Lxyz=data[index][14]
			break;
		case 9:
			return data[ind][15]; //Lxzz=data[index][15]
			break;
		case 8:
			return data[ind][16]; //Lyyy=data[index][16]
			break;
		case 12:
			return data[ind][17]; //Lyyz=data[index][17]
			break;
		case 18:
			return data[ind][18]; //Lyzz=data[index][18]
			break;
		case 27:
			return data[ind][19]; //Lzzz=data[index][19]
			break;
		}

		// should not get past switch statement
		assert(0);
		return NaN;
	}

	//---------------------------------------------------
	// increment tensor indices for any symmetric tensor
	//---------------------------------------------------
	static void increment_tensor_indices(unsigned tensor_indices[], const unsigned rank) {
		PRINT_ASSERT(rank,>=,0);
		PRINT_ASSERT(rank,<,nranks);
		if (rank == 0)
			return;
		tensor_indices[0]++;
		for (int which_index = 0; which_index < rank - 1; which_index++) {
			if (tensor_indices[which_index] > index_range - 1) {
				tensor_indices[which_index + 1]++;
				for (int which_index_2 = 0; which_index_2 < which_index + 1;
						which_index_2++)
					tensor_indices[which_index_2] = tensor_indices[which_index + 1];
			}
		}
	}

	//--------------------------------------------------------------
	// Initialization and Allocation
	//--------------------------------------------------------------

	void init(const vector<Axis>& spatial_axes, const Axis& nu_grid) {
		PRINT_ASSERT(spatial_axes.size(),==,ndims_spatial);
		vector<Axis> axes(ndims_spatial+1);

		// spatial axes
		for(int i=0; i<ndims_spatial; i++) axes[i] = spatial_axes[i];

		// frequency axis
		axes[nuGridIndex] = nu_grid;

		// initialize the moments
		data.set_axes(axes);
		wipe();
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
	void count(const EinsteinHelper* eh, const double E) {
		PRINT_ASSERT(E, >=, 0);
		PRINT_ASSERT(E, <, INFINITY);
		double D[3] = {eh->kup_tet[0], eh->kup_tet[1], eh->kup_tet[2]};
		Metric::normalize_Minkowski<3>(D);
		
		unsigned data_indices[data.Ndims()];
		for(int i=0; i<ndims_spatial; i++) data_indices[i] = eh->dir_ind[i];
		data_indices[nuGridIndex] = eh->dir_ind[NDIMS];

		// increment moments
		Tuple<double, n_total_elements> tmp;
		unsigned tuple_index=0;
		for(unsigned rank = 0; rank<nranks; rank++) {
			unsigned tensor_indices[rank];
			for(unsigned r = 0; r<rank; r++) tensor_indices[r] = 0;
			for(unsigned i=0; i<n_independent_elements[rank]; i++) {
				tmp[tuple_index] = E;
				for(int r = 0; r<rank; r++) tmp[tuple_index] *= D[tensor_indices[r]];
				increment_tensor_indices(tensor_indices, rank);
				tuple_index++;
			}
		}
		data.add(data_indices, tmp);
	}

	double reconstruct_f(const unsigned dir_ind[ndims_spatial+1], const double k[3]) const{
		unsigned index = data.direct_index(dir_ind);
		Tuple<double,n_total_elements> M = data[index];
		const double x=k[0], y=k[1], z=k[2];
		const double E=M[0];
		const double Fx=M[1], Fy=M[2], Fz=M[3];
		const double Pxx=M[4], Pxy=M[5], Pxz=M[6], Pyy=M[7], Pyz=M[8], Pzz=M[9];
		const double Lxxx=M[10], Lxxy=M[11], Lxxz=M[12], Lxyy=M[13], Lxyz=M[14];
		const double Lxzz=M[15], Lyyy=M[16], Lyyz=M[17], Lyzz=M[18], Lzzz=M[19];

		double f = E;
		f += 3. * (x*Fx + y*Fy * z*Fz);
		f += 15./2. * (x*x*Pxx + y*y*Pyy + z*z*Pzz
				+ 2.*(x*y*Pxy + x*z*Pxz + y*z*Pyz) - E/3.);
		f += 35./2. * (x*x*x*Lxxx + y*y*y*Lyyy + z*z*z*Lzzz
				+ 3. * (x*x*(z*Lxxz + y*Lxxy) + y*y*(x*Lxyy + z*Lyyz) + z*z*(x*Lxzz + y*Lyzz))
				+ 6.*x*y*z*Lxyz - 3./5.*(x*Fx + y*Fy + z*Fz));

		return f/(4.*physical_constants::pi);
	}


	void rescale(const double r) {
		for(unsigned i=0; i<data.size(); i++) data[i] *= r;
	}
	void rescale_spatial_point(const unsigned dir_ind[ndims_spatial], const double r){
		unsigned all_indices[ndims_spatial+1];
		for(unsigned i=0; i<ndims_spatial; i++) all_indices[i] = dir_ind[i];
		all_indices[ndims_spatial  ] = 0;
		unsigned base_ind = data.direct_index(all_indices);
		unsigned nbins = data.axes[ndims_spatial].size();
		for(unsigned i=0; i<nbins; i++){
			data.y0[base_ind+i] *= r;
		}
	}

	//--------------------------------------------------------------
	// MPI scatter the spectrum contents.
	// Must rescale zone stop list to account for number of groups
	//--------------------------------------------------------------
	void mpi_sum_scatter(vector<unsigned>& zone_stop_list){
		unsigned ngroups = data.axes[nuGridIndex].size();
		vector<unsigned> stop_list = zone_stop_list;
		for(unsigned i=0; i<stop_list.size(); i++) stop_list[i] *= ngroups;
		data.mpi_sum_scatter(stop_list);
	}
	void mpi_sum(){
		data.mpi_sum();
	}

	//--------------------------------------------------------------
	// Write data to specified location in an HDF5 file
	//--------------------------------------------------------------
	void write_hdf5_data(H5::H5File file,const string name) {
		data.write_HDF5(file, name);
	}

	//--------------------------------------------------------------
	// Write distribution function coordinates to an HDF5 file
	//--------------------------------------------------------------
	void write_hdf5_coordinates(H5::H5File file,
			const string name) const {
		// useful quantities
		hsize_t dims[1];
		H5::DataSpace dataspace;
		H5::DataSet dataset;
		H5::Group group;
		vector<float> tmp;

		// set up the dfunc group
		std::stringstream datasetname, indicesname;

		// SET UP DATASPACE FOR EACH MOMENT
		for (unsigned rank=0; rank<nranks; rank++) {
			// prep the filenames
			indicesname.str("");
			indicesname << name;
			indicesname << "_moment_" << rank << "_indices";

			// set up the database with the indices
			hsize_t indices_dimensions[] = { n_independent_elements[rank], rank };
			dataspace = H5::DataSpace(2, indices_dimensions);
			dataset = file.createDataSet(indicesname.str(), H5::PredType::STD_I32LE,
					dataspace);
			unsigned tensor_indices[rank];
			for (int r = 0; r < rank; r++)
				tensor_indices[r] = 0;
			int indices2D[n_independent_elements[rank]][rank];
			for (int i = 0; i < n_independent_elements[rank]; i++) {
				for (int r = 0; r < rank; r++)
					indices2D[i][r] = tensor_indices[r];
				increment_tensor_indices(tensor_indices, rank);
			}
			dataset.write(indices2D, H5::PredType::STD_I32LE);
			dataset.close();

		}
	}
	void add_isotropic(const unsigned dir_ind[NDIMS+1], const double E){
		unsigned indices[data.Ndims()];
		for(int i=0; i<ndims_spatial; i++) indices[i] = dir_ind[i];
		indices[nuGridIndex] = dir_ind[ndims_spatial];

		// increment moments
		Tuple<double, n_total_elements> tmp;
		tmp = 0;
		tmp[0] = E;
		tmp[4] = E/3.; // xx
		tmp[7] = E/3.; // yy
		tmp[9] = E/3.; // zz
		data.add(indices, tmp);
	}

	double total() const{
		double result=0;
		for(unsigned i=0; i<data.size(); i++)
			result += data[i][0];
		return result;
	}

	void annihilation_rate(
			const unsigned dir_ind[NDIMS],       // spatial directional indices for the zone we're getting the rate at
			const SpectrumArray* in_dist,  // erg/ccm (integrated over angular bin and energy bin)
			const vector< vector<vector<double> > >& phi, // cm^3/s [order][igin][igout]
			const unsigned weight,
			Tuple<double,4>& fourforce) const{

		const MomentSpectrumArray<NDIMS>* nubar_dist = (MomentSpectrumArray<NDIMS>*)in_dist;
		PRINT_ASSERT(phi.size(),>=,2);

		unsigned tmp_ind[NDIMS+1];
		for(unsigned i=0; i<NDIMS; i++) tmp_ind[i] = dir_ind[i];
		tmp_ind[NDIMS] = 0;
		const unsigned base_ind = direct_index(tmp_ind);
		const Axis* nu_axis = &(data.axes[nuGridIndex]);
		const unsigned nnu = nu_axis->size();

		for(unsigned i=0; i<nnu; i++){
			double avg_e = pc::h * nu_axis->mid[i];
			for(unsigned j=0; j<nnu; j++){
				double avg_ebar = pc::h * nu_axis->mid[j];
				double eebar = avg_e*avg_ebar;
				double phi2 = -1./5. * (phi[0][i][j] + 3.*phi[1][i][j]);

				// basic spherically symmetric
				double tmp0 = getE(i+base_ind)*nubar_dist->getE(j+base_ind) / eebar; // #/cm^6
				PRINT_ASSERT(tmp0,>=,0);
				double tmp1=0, tmp2=0;
				for(unsigned k=0; k<3; k++){
					tmp1 += getF(i+base_ind,k)*nubar_dist->getF(j+base_ind,k) / eebar; // #/cm^6
					for(unsigned l=0; l<3; l++) tmp2 += 3.*getP(i+base_ind,k,l)*nubar_dist->getP(j+base_ind,k,l)/eebar;
					tmp2 -= tmp0;
					tmp2 *= 0.5;
					PRINT_ASSERT(abs(getF(i+base_ind,k)),<=,getE(i+base_ind));
					PRINT_ASSERT(abs(nubar_dist->getF(i+base_ind,k)),<=,nubar_dist->getE(i+base_ind));
				}
				double this_dep = (avg_e + avg_ebar) * (0.5*phi[0][i][j]*tmp0 + 1.5*phi[1][i][j]*tmp1 + 2.5*phi2*tmp2); // erg/cm^3/s
				fourforce[3] += this_dep;

				// space components
				for(unsigned a=0; a<3; a++){
					tmp0 = avg_e    * getF(i+base_ind,a) * nubar_dist->getE(j+base_ind)   / eebar + // erg/cm^6
						   avg_ebar * getE(i+base_ind)   * nubar_dist->getF(j+base_ind,a) / eebar;
					tmp1=0;
					tmp2=0;
					for(unsigned b=0; b<3; b++){
						tmp1 += avg_e    * getP(i+base_ind,a,b) * nubar_dist->getF(j+base_ind,  b) / eebar; // erg/cm^6
						tmp1 += avg_ebar * getF(i+base_ind,  b) * nubar_dist->getP(j+base_ind,a,b) / eebar;
						for(unsigned c=0; c<3; c++){
							tmp2 += avg_e    * 3.*getL(i+base_ind,a,b,c)*nubar_dist->getP(j+base_ind,  b,c)/eebar;
							tmp2 += avg_ebar * 3.*getP(i+base_ind,  b,c)*nubar_dist->getL(j+base_ind,a,b,c)/eebar;
						}
						tmp2 -= tmp0;
						tmp2 *= 0.5;
					}
					fourforce[a] += 0.5*phi[0][i][j]*tmp0 + 1.5*phi[1][i][j]*tmp1 + 2.5*phi2*tmp2; // erg/ccm/s
				}
			}
		}

		// sanity checks
		for(unsigned i=0; i<4; i++) fourforce[i] /= weight;
		//PRINT_ASSERT(fourforce[3],>=,0); // not always true, since moment reconstruction and annihilation kernel are not exact
	}

};

#endif
