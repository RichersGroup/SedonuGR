#include <limits>
#include <omp.h>
#include <mpi.h>
#include <iostream>
#include <math.h>
#include <string.h>
#include <vector>
#include <fstream>
#include "global_options.h"
#include "spectrum_array.h"

//--------------------------------------------------------------
// Initialization and Allocation
//--------------------------------------------------------------
void spectrum_array::init(const std::vector<double> w,
		const int n_mu, const int n_phi)
{
	// assign wave grid
	double w_start = w[0];
	double w_stop  = w[1];
	double w_del   = w[2];
	nu_grid.init(w_start,w_stop,w_del);
	int n_wave   = nu_grid.size();

	// asign mu grid
	mu_grid.init(-1,1,n_mu);

	// asign phi grid
	phi_grid.init(-pc::pi,pc::pi,n_phi);

	// index parameters
	unsigned n_elements  = n_wave*n_mu*n_phi;

	// allocate memory
	flux.resize(n_elements);

	// clear
	wipe();
}


//--------------------------------------------------------------
// Initialization and Allocation
//--------------------------------------------------------------
void spectrum_array::init(const locate_array wg,
		const locate_array mg, const locate_array pg)
{
	// initialize locate arrays by swapping with the inputs
	nu_grid.copy(wg);
	mu_grid.copy(mg);
	phi_grid.copy(pg);

	int n_wave   = nu_grid.size();
	int n_mu     = mu_grid.size();
	int n_phi    = phi_grid.size();

	// index parameters
	unsigned n_elements  = n_wave*n_mu*n_phi;

	// allocate memory
	flux.resize(n_elements);

	// clear
	wipe();
}


//--------------------------------------------------------------
// Functional procedure: Wipe
//--------------------------------------------------------------
void spectrum_array::wipe()
{
	for(unsigned i=0;i<flux.size();i++) flux[i] = 0;
}



//--------------------------------------------------------------
// handles the indexing: should be called in this order
//    group, mu, phi
//--------------------------------------------------------------
unsigned spectrum_array::index(const unsigned nu_bin, const unsigned mu_bin, const unsigned phi_bin) const
{
	unsigned n_phi = phi_grid.size();
	unsigned n_nu  =  nu_grid.size();
	unsigned n_mu  =  mu_grid.size();
	assert(nu_bin < n_nu);
	assert(mu_bin < n_mu);
	assert(phi_bin < n_phi);
	int a3 = n_phi;
	int a2 = n_mu*a3;
	const unsigned ind = nu_bin*a2 + mu_bin*a3 + phi_bin;
	assert(ind < n_phi*n_nu*n_mu);
	return ind;
}

//----------------
// get bin indices
//----------------
unsigned spectrum_array::nu_bin(const unsigned index) const{
	assert(index<flux.size());
	unsigned result = index / (phi_grid.size()*mu_grid.size());
	assert(result<nu_grid.size());
	return result;
}
unsigned spectrum_array::mu_bin(const unsigned index) const{
	assert(index<flux.size());
	double result =  (index%(phi_grid.size()*mu_grid.size())) / phi_grid.size();
	assert(result<mu_grid.size());
	return result;
}
unsigned spectrum_array::phi_bin(const unsigned index) const{
	assert(index<flux.size());
	double result = index % phi_grid.size();
	assert(result<phi_grid.size());
	return result;
}

//--------------------
// get centers of bins
//--------------------
double spectrum_array::nu_center(const unsigned index) const{
	assert(index<flux.size());
	return nu_grid.center(nu_bin(index));
}
double spectrum_array::mu_center(const unsigned index) const{
	assert(index<flux.size());
	return mu_grid.center(mu_bin(index));
}
double spectrum_array::phi_center(const unsigned index) const{
	assert(index<flux.size());
	return phi_grid.center(phi_bin(index));
}

//-------
// getter
//-------
double spectrum_array::get(const int index) const{
	return flux[index];
}

//----------------------
// get size of spectrum
//----------------------
unsigned spectrum_array::size() const{
	return flux.size();
}

//--------------------------------------------------------------
// count a particle
////--------------------------------------------------------------
void spectrum_array::count(const vector<double>& D, const double nu, const double E)
{
	assert(D.size()==3);
	assert(E>=0);
	double mu  = D[2];           // component along z axis
	mu = max(-1.0,mu);
	mu = min( 1.0,mu);
	double phi = atan2(D[1],D[0]);  // projection into x-y plane
	if(phi< -pc::pi) phi += 2.0*pc::pi;
	if(phi>  pc::pi) phi -= 2.0*pc::pi;

	// if off the LEFT of mu/phi grids, just return without counting
	if ((mu<mu_grid.min) || (phi<phi_grid.min)){
		cout << "Lost a particle off the left of the spectrum array!" << endl;
		return;
	}

	// locate bin number in all dimensions.
	unsigned nu_bin  = nu_grid.locate(nu);
	unsigned mu_bin  =   mu_grid.locate(mu);
	unsigned phi_bin =  phi_grid.locate(phi);

	// if off the RIGHT of mu/phi grids, just return without counting
	if((mu_bin ==   mu_grid.size()) || (phi_bin ==  phi_grid.size())){
		cout << "Lost a particle off the right of the spectrum array!" << endl;
		return;
	}

	// if off RIGHT of wavelength grid, store in last bin (LEFT is accounted for by locate)
	if (nu_bin == nu_grid.size()) nu_bin--;

	// add to counters
	int ind = index(nu_bin,mu_bin,phi_bin);

    #pragma omp atomic
	flux[ind]  += E;
}

double spectrum_array::integrate() const{
	double integral = 0;
	for(unsigned nu_bin=0; nu_bin<nu_grid.size(); nu_bin++){
		for(unsigned mu_bin=0; mu_bin<mu_grid.size(); mu_bin++){
			for(unsigned phi_bin=0; phi_bin<phi_grid.size(); phi_bin++){
				int ind = index(nu_bin,mu_bin,phi_bin);
				integral += flux[ind];
			}
		}
	}
	return integral;
}

// integrate over direction, but also divide by bin width
void spectrum_array::integrate_over_direction(vector<double>& integral) const{
	integral = vector<double>(nu_grid.size(),0);
	for(unsigned nu_bin=0; nu_bin<nu_grid.size(); nu_bin++){
		for(unsigned mu_bin=0; mu_bin<mu_grid.size(); mu_bin++){
			for(unsigned phi_bin=0; phi_bin<phi_grid.size(); phi_bin++){
				int ind = index(nu_bin,mu_bin,phi_bin);
				integral[nu_bin] += flux[ind]/(nu_grid.delta(nu_bin)*pc::h_MeV);
			}
		}
	}
}

//--------------------------------------------------------------
// print out
//--------------------------------------------------------------
void spectrum_array::print(const int iw, const int species) const
{
	ofstream outf;
	string filename = "spectrum_species"+to_string(species);
	transport::open_file(filename.c_str(),iw,outf);

	unsigned n_nu  =  nu_grid.size();
	unsigned n_mu  =  mu_grid.size();
	unsigned n_phi = phi_grid.size();

	outf << "# " << "n_nu:" << n_nu << " " << "n_mu:" << n_mu << " " << "n_phi:" << n_phi << endl;
	outf << "# ";
	outf << (n_nu >1 ? "frequency(Hz) " : "");
	outf << (n_mu >1 ? "mu "            : "");
	outf << (n_phi>1 ? "phi "           : "");
	outf << "integrated_flux(erg) counts" << endl;

	for (unsigned k=0;k<n_mu;k++)
		for (unsigned m=0;m<n_phi;m++)
			for (unsigned j=0;j<n_nu;j++)
			{
				int id = index(j,k,m);
				if (n_nu > 1) outf <<  nu_grid.center(j) << " ";
				if (n_mu > 1) outf <<  mu_grid.center(k) << " ";
				if (n_phi> 1) outf << phi_grid.center(m) << " ";

				// the delta is infinity if the bin is a catch-all.
				double wdel = nu_grid.delta(j);
				double norm = n_mu*n_phi
						* ( wdel < numeric_limits<double>::infinity() ? wdel : 1 );
				outf << flux[id]/norm << endl;
			}
	outf.close();
}


void  spectrum_array::rescale(double r)
{
	for(unsigned i=0;i<flux.size();i++) flux[i] *= r;
}


//--------------------------------------------------------------
// MPI average the spectrum contents
//--------------------------------------------------------------
// only process 0 gets the reduced spectrum to print
void spectrum_array::MPI_average()
{
	const unsigned n_elements = nu_grid.size()*mu_grid.size()*phi_grid.size();

	// average the flux (receive goes out of scope after section)
	vector<double> receive;
	receive.resize(n_elements);
	MPI_Allreduce(&flux.front(), &receive.front(), n_elements, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	for(unsigned i=0; i<flux.size(); i++) flux[i] = receive[i];//flux.swap(receive);

	// only have the receiving ID do the division
	int myID, mpi_procs;
	MPI_Comm_size( MPI_COMM_WORLD, &mpi_procs );
	MPI_Comm_rank( MPI_COMM_WORLD, &myID      );
	rescale(1./(double)mpi_procs);
}
