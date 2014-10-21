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
	nu_grid.swap(wg);
	mu_grid.swap(mg);
	phi_grid.swap(pg);

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
int spectrum_array::index(const int nu_bin, const int mu_bin, const int phi_bin) const
{
	assert(nu_bin>=0);
	assert(mu_bin>=0);
	assert(phi_bin>=0);
	int n_phi = phi_grid.size();
	int n_nu  =  nu_grid.size();
	int n_mu  =  mu_grid.size();
	int a3 = n_phi;
	int a2 = n_mu*a3;
	const int ind = nu_bin*a2 + mu_bin*a3 + phi_bin;
	assert(ind >= 0);
	assert(ind < n_phi*n_nu*n_mu);
	return ind;
}


//--------------------------------------------------------------
// count a particle
////--------------------------------------------------------------
void spectrum_array::count(const particle* p, const double E)
{
	assert(p->D.size()==3);
	assert(E>=0);
	double mu  = p->D[2];           // component along z axis
	double phi = atan2(p->D[1],p->D[0]);  // projection into x-y plane

	// if off the LEFT of mu/phi grids, just return without counting
	if ((mu<mu_grid.min) || (phi<phi_grid.min)) return;

	// locate bin number in all dimensions.
	unsigned nu_bin  = nu_grid.locate(p->nu);
	unsigned mu_bin  =   mu_grid.locate(mu);
	unsigned phi_bin =  phi_grid.locate(phi);

	// if off the RIGHT of mu/phi grids, just return without counting
	if((mu_bin ==   mu_grid.size()) || (phi_bin ==  phi_grid.size())) return;

	// if off RIGHT of wavelength grid, store in last bin (LEFT is accounted for by locate)
	if (nu_bin == nu_grid.size()) nu_bin--;

	// add to counters
	int ind = index(nu_bin,mu_bin,phi_bin);

    #pragma omp atomic
	flux[ind]  += E;
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
	flux.swap(receive);

	// only have the receiving ID do the division
	int myID, mpi_procs;
	MPI_Comm_size( MPI_COMM_WORLD, &mpi_procs );
	MPI_Comm_rank( MPI_COMM_WORLD, &myID      );
	if(myID == 0) rescale(1./(double)mpi_procs);
}
