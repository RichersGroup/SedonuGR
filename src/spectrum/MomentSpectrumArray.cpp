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
#include "MomentSpectrumArray.h"
#include "Transport.h"

using namespace std;
namespace pc = physical_constants;

//---------------------------------------------------
// increment tensor indices for any symmetric tensor
//---------------------------------------------------
void increment_tensor_indices(int tensor_indices[], int rank, int length) {
	PRINT_ASSERT(rank, >=, 0);
	PRINT_ASSERT(length, >, 0);
	if (rank == 0)
		return;
	tensor_indices[0]++;
	for (int which_index = 0; which_index < rank - 1; which_index++) {
		if (tensor_indices[which_index] > length - 1) {
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
void MomentSpectrumArray::init(const LocateArray wg, const int max_rank) {
	// initialize locate arrays by swapping with the inputs
	nu_grid.copy(wg);

	// determine the number of independent elements
	const int length = 3;
	int n_independent_elements[max_rank + 1];
	n_independent_elements[0] = 1;
	for (int rank = 1; rank < max_rank + 1; rank++) {
		n_independent_elements[rank] = 0;
		int tensor_indices[rank];
		for (int i = 0; i < rank; i++)
			tensor_indices[i] = 0;
		while (tensor_indices[rank - 1] < length) {
			n_independent_elements[rank]++;
			increment_tensor_indices(tensor_indices, rank, length);
		}
	}

	// initialize the moments
	moments.resize(nu_grid.size());
	for (int inu = 0; inu < nu_grid.size(); inu++) {
		moments[inu].resize(max_rank + 1);
		for (int rank = 0; rank < max_rank + 1; rank++) {
			moments[inu][rank].resize(n_independent_elements[rank]);
		}
	}

	// clear the data
	wipe();
}

//--------------------------------------------------------------
// Functional procedure: Wipe
//--------------------------------------------------------------
void MomentSpectrumArray::wipe() {
	for (int group = 0; group < moments.size(); group++)
		for (int rank = 0; rank < moments[group].size(); rank++)
			for (int i = 0; i < moments[group][rank].size(); i++)
				moments[group][rank][i] = 0;
}

//----------------------
// get bin central value
//----------------------
double MomentSpectrumArray::nu_bin_center(const unsigned index) const {
	PRINT_ASSERT(index, <, nu_grid.size());
	return nu_grid.center(index);
}

//----------------------
// get size of spectrum
//----------------------
unsigned MomentSpectrumArray::n_groups() const {
	return moments.size();
}
unsigned MomentSpectrumArray::n_ranks() const {
	return moments[0].size();
}
unsigned MomentSpectrumArray::n_elements(const int rank) const {
	PRINT_ASSERT(rank, <, moments[0][rank].size());
	return moments[0][rank].size();
}

//--------------------------------------------------------------
// count a particle
////--------------------------------------------------------------
void MomentSpectrumArray::count(const double D[3], const int Dsize,
		const double nu, const double E) {
	PRINT_ASSERT(Dsize, ==, 3);
	PRINT_ASSERT(E, >=, 0);
	PRINT_ASSERT(E, <, INFINITY);

	// locate bin number
	unsigned nu_bin = nu_grid.locate(nu);
	if (nu_bin == nu_grid.size()) nu_bin--;

	// increment moments
	double tmp;
	for (int rank = 0; rank<n_ranks(); rank++) {
		int tensor_indices[rank];
		for (int r = 0; r<rank; r++) tensor_indices[r] = 0;
		
		for (int i = 0; i<n_elements(rank); i++) {
			tmp = E;
			for (int r = 0; r<rank; r++) tmp *= D[tensor_indices[r]];
			increment_tensor_indices(tensor_indices, rank, 3);

			#pragma omp atomic
			moments[nu_bin][rank][i] += tmp;
		}
	}
}

double MomentSpectrumArray::average_nu() const {
	double integral_E = 0;
	double integral_N = 0;
	for (unsigned nu_bin = 0; nu_bin<n_groups(); nu_bin++) {
		integral_E += moments[nu_bin][0][0];
		integral_N += moments[nu_bin][0][0] / nu_grid.center(nu_bin);
	}
	double result = integral_E / integral_N;
	if (result >= 0)
		return result;
	else
		return 0;
}

double MomentSpectrumArray::integrate() const {
	double integral = 0;
	for (unsigned nu_bin = 0; nu_bin<n_groups(); nu_bin++)
		integral += moments[nu_bin][0][0];
	return integral;
}

// integrate over direction
void MomentSpectrumArray::integrate_over_direction(
		vector<double>& integral) const {
	integral = vector<double>(nu_grid.size(), 0);
	for (unsigned nu_bin = 0; nu_bin<n_groups(); nu_bin++) {
		integral[nu_bin] = moments[nu_bin][0][0];
	}
}

//--------------------------------------------------------------
// print out
//--------------------------------------------------------------
void MomentSpectrumArray::print(const int iw, const int species) const {
	ofstream outf;
	stringstream speciesstream;
	speciesstream << species;
	string prefix = "spectrum_species" + speciesstream.str();
	string filename = Transport::filename(prefix.c_str(), iw, ".dat");
	outf.open(filename.c_str());

	// write the header
	outf << "# " << "n_nu:" << n_groups() << " n_ranks:" << n_ranks() << endl;
	outf << "# ";
	outf << "frequency(Hz) delta(Hz) M0(erg/s/Hz)";
	for (int rank = 1; rank < n_ranks(); rank++) {
		int tensor_indices[rank];
		for (int r = 0; r < rank; r++)
			tensor_indices[r] = 0;
		while (tensor_indices[rank - 1] < 3) {
			outf << " M" << rank << ":";
			for (int r = 0; r < rank; r++)
				outf << tensor_indices[r];
		}
	}
	outf << endl;

	// write the data
	for (unsigned group = 0; group < n_groups(); group++) {
		outf << nu_grid.center(group) << " " << nu_grid.delta(group) << " ";
		for (unsigned rank = 0; rank < n_ranks(); rank++)
			for (unsigned i = 0; i < n_elements(rank); i++) {
				// the delta is infinity if the bin is a catch-all.
				double wdel = (
						nu_grid.delta(group)
								< numeric_limits<double>::infinity() ?
								nu_grid.delta(group) : 1);
				outf << moments[group][rank][i] / wdel << endl;
			}
	}
	outf.close();
}

void MomentSpectrumArray::rescale(double r) {
	for (int group = 0; group < n_groups(); group++)
		for (int rank = 0; rank < n_ranks(); rank++)
			for (int i = 0; i < n_elements(rank); i++)
				moments[group][rank][i] *= r;
}

//--------------------------------------------------------------
// MPI average the spectrum contents
//--------------------------------------------------------------
// only process 0 gets the reduced spectrum to print
void MomentSpectrumArray::MPI_average() {
	vector<double> receive;

	for (int group = 0; group < n_groups(); group++)
		for (int rank = 0; rank < n_ranks(); rank++) {
			receive.resize(n_elements(rank));
			MPI_Allreduce(&moments[group][rank].front(), &receive.front(),
					n_elements(rank), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			for (unsigned i = 0; i < n_elements(rank); i++)
				moments[group][rank][i] = receive[i]; //flux.swap(receive);
		}

	// only have the receiving ID do the division
	int myID, mpi_procs;
	MPI_Comm_size( MPI_COMM_WORLD, &mpi_procs);
	MPI_Comm_rank( MPI_COMM_WORLD, &myID);
	rescale(1. / (double) mpi_procs);
}

//--------------------------------------------------------------
// Write data to specified location in an HDF5 file
//--------------------------------------------------------------
void MomentSpectrumArray::write_hdf5_data(H5::H5File file, const int s,
		const int dir_ind[], const hsize_t n_spatial_dims) const {

	for (int rank = 0; rank < n_ranks(); rank++) {
		// prep the filenames
		stringstream datasetname, indicesname;
		string dfunc_dir = "distribution(erg|ccm,lab)";
		datasetname << dfunc_dir << "/M" << rank;
		vector<float> tmp(n_elements(rank), -1.0);

		// get the dataset
		H5::DataSet dataset = file.openDataSet(datasetname.str());

		// get the dataspace of the dataset
		H5::DataSpace dataspace = dataset.getSpace();
		hsize_t n_total_dims = dataspace.getSimpleExtentNdims();
		hsize_t dfunc_dims[n_spatial_dims + 3];
		hsize_t start[n_spatial_dims + 3];
		dataspace.getSelectBounds(start, dfunc_dims);
		for (int i = 0; i < n_total_dims; i++)
			dfunc_dims[i] += 1;
		PRINT_ASSERT(dataspace.getSimpleExtentNdims(), ==, n_spatial_dims + 3);
		PRINT_ASSERT(dfunc_dims[n_spatial_dims + 1], ==, n_groups());
		PRINT_ASSERT(dfunc_dims[n_spatial_dims + 2], ==, n_elements(rank));

		for (int group = 0; group < n_groups(); group++) {
			// set up subspace offset
			hsize_t offset[n_total_dims];
			for (unsigned i = 0; i < n_spatial_dims; i++)
				offset[i] = dir_ind[i];
			offset[n_spatial_dims] = s;
			offset[n_spatial_dims + 1] = group;
			offset[n_spatial_dims + 2] = 0;

			// set up subspace stride
			hsize_t stride[n_total_dims];
			for (unsigned i = 0; i < n_total_dims; i++)
				stride[i] = 1;

			// set up subspace block
			hsize_t block[n_total_dims];
			for (unsigned i = 0; i < n_total_dims; i++)
				block[i] = 1;

			// set up local spectrum dimensions
			hsize_t spectrum_dims[n_total_dims];
			for (unsigned i = 0; i < n_spatial_dims; i++)
				spectrum_dims[i] = 1;
			spectrum_dims[n_spatial_dims] = 1;
			spectrum_dims[n_spatial_dims + 1] = 1;
			spectrum_dims[n_spatial_dims + 2] = n_elements(rank);

			// set dataspace
			dataspace.selectHyperslab(H5S_SELECT_SET, &spectrum_dims[0],
					&offset[0], &stride[0], &block[0]);

			// define the memory dataspace
			hsize_t mem_length = n_elements(rank);
			H5::DataSpace memspace(1, &mem_length, NULL);

			// write the data (converting to single precision)
			// assumes phi increases fastest, then mu, then nu
			for (unsigned i = 0; i < n_elements(rank); i++)
				tmp[i] = moments[group][rank][i];
			dataset.write(&tmp[0], H5::PredType::IEEE_F32LE, memspace,
					dataspace);
		}
		dataset.close();
	}
}

//--------------------------------------------------------------
// Write distribution function coordinates to an HDF5 file
//--------------------------------------------------------------
void MomentSpectrumArray::write_hdf5_coordinates(H5::H5File file,
		const Grid* grid) const {
	// useful quantities
	hsize_t dims[1];
	H5::DataSpace dataspace;
	H5::DataSet dataset;
	H5::Group group;
	vector<float> tmp;

	// set up the dfunc group
	string dfunc_dir = "distribution(erg|ccm,lab)";
	group = file.createGroup(dfunc_dir);
	stringstream datasetname, indicesname;

	// write nu_grid
	tmp = vector<float>(nu_grid.size() + 1, 0.0);
	dims[0] = nu_grid.size() + 1;
	dataspace = H5::DataSpace(1, dims);

	// write the frequency grid
	datasetname << dfunc_dir << "/distribution_frequency_grid(Hz,lab)";
	dataset = file.createDataSet(datasetname.str(), H5::PredType::IEEE_F32LE,
			dataspace);
	tmp[0] = nu_grid.min;
	for (unsigned i = 1; i < nu_grid.size() + 1; i++)
		tmp[i] = nu_grid[i - 1];
	dataset.write(&tmp[0], H5::PredType::IEEE_F32LE);
	dataset.close();

	// set up list of dimensions in each direction
	hsize_t zdims[grid->dimensionality()];
	grid->dims(zdims, grid->dimensionality());
	vector<hsize_t> dims_plus3(grid->dimensionality() + 3, 0);
	for (unsigned i = 0; i < grid->dimensionality(); i++)
		dims_plus3[i] = zdims[i]; // number of spatial bins
	dims_plus3[grid->dimensionality()] = grid->z[0].distribution.size(); // number of species
	dims_plus3[grid->dimensionality() + 1] = n_groups(); // number of energy bins
	dims_plus3[grid->dimensionality() + 2] = 0; // number of elements in moment

	// SET UP DATASPACE FOR EACH MOMENT
	for (unsigned rank = 0; rank < n_ranks(); rank++) {
		// prep the filenames
		datasetname.str("");
		datasetname << dfunc_dir << "/M" << rank;
		indicesname.str("");
		indicesname << dfunc_dir << "/M" << rank << "indices";

		// set the HDF5 dataspace
		dims_plus3[grid->dimensionality() + 2] = n_elements(rank);
		dataspace = H5::DataSpace(grid->dimensionality() + 3, &dims_plus3[0]);
		dataset = file.createDataSet(datasetname.str(),
				H5::PredType::IEEE_F32LE, dataspace);

		// set up the database with the indices
		hsize_t indices_dimensions[] = { n_elements(rank), rank };
		dataspace = H5::DataSpace(2, indices_dimensions);
		dataset = file.createDataSet(indicesname.str(), H5::PredType::STD_I32LE,
				dataspace);
		int tensor_indices[rank];
		for (int r = 0; r < rank; r++)
			tensor_indices[r] = 0;
		int indices2D[n_elements(rank)][rank];
		for (int i = 0; i < n_elements(rank); i++) {
			for (int r = 0; r < rank; r++)
				indices2D[i][r] = tensor_indices[r];
			increment_tensor_indices(tensor_indices, rank, 3);
		}
		dataset.write(indices2D, H5::PredType::STD_I32LE);
		dataset.close();

	}
}

void MomentSpectrumArray::write_header(ofstream& outf) const {
	outf << "frequency(Hz) delta(Hz) M0(erg/s/Hz/sr)";
	for (int rank = 1; rank < n_ranks(); rank++) {
		int tensor_indices[rank];
		for (int r = 0; r < rank; r++)
			tensor_indices[r] = 0;
		while (tensor_indices[rank - 1] < 3) {
			outf << " M" << rank << ":";
			for (int r = 0; r < rank; r++)
				outf << tensor_indices[r];
		}
	}

	for (int group = 0; group < n_groups(); group++)
		for (int rank = 0; rank < n_ranks(); rank++) {
			if (rank > 0) {
				outf << "i";
				int tensor_indices[rank];
				for (int r = 0; r < rank; r++)
					tensor_indices[r] = 0;
				while (tensor_indices[rank - 1] < 3) {
					outf << "g" << group << "M" << rank << ":";
					for (int r = 0; r < rank; r++)
						outf << tensor_indices[r];
				}
			} else
				outf << "g" << group << "M0";
			outf << "edens(erg/ccm)  ";
		}
}

void MomentSpectrumArray::write_line(ofstream& outf) const {
	for (int group = 0; group < moments.size(); group++)
		for (int rank = 0; rank < moments[group].size(); rank++)
			for (int i = 0; i < moments[group][rank].size(); i++)
				outf << moments[group][rank][i] << "\t";
}
