#pragma warning disable 161
#include <limits>
#include <omp.h>
#include <mpi.h>
#include <stdio.h>
#include <iostream>
#include <math.h>
#include <string.h>
#include <vector>
#include "spectrum_array.h"
#include "physical_constants.h"

namespace pc = physical_constants;

//--------------------------------------------------------------
// Constructors
//--------------------------------------------------------------
spectrum_array::spectrum_array()
{
  strcpy(name,DEFAULT_NAME);
  a1=0; a2=0; a3=0;
  n_elements=0;
}

void spectrum_array::set_name(const char *n)
{
  strcpy(name,n);
}


//--------------------------------------------------------------
// Initialization and Allocation
//--------------------------------------------------------------
void spectrum_array::init(std::vector<double> t, std::vector<double> w,
			  int n_mu, int n_phi)
{
  // assign time grid
  double t_start = t[0];
  double t_stop  = t[1];
  double t_del   = t[2];
  time_grid.init(t_start,t_stop,t_del);
  int n_times  = time_grid.size();

  // assign wave grid
  double w_start = w[0];
  double w_stop  = w[1];
  double w_del   = w[2];
  wave_grid.init(w_start,w_stop,w_del);
  int n_wave   = wave_grid.size();

  // asign mu grid
  mu_grid.init(-1,1,n_mu);

  // asign phi grid
  phi_grid.init(-pc::pi,pc::pi,n_phi);

  // index parameters
  n_elements  = n_times*n_wave*n_mu*n_phi;
  a3 = n_phi;
  a2 = n_mu*a3;
  a1 = n_wave*a2;

  // allocate memory
  click.resize(n_elements);
  flux.resize(n_elements);

  // clear 
  wipe();
}


//--------------------------------------------------------------
// Initialization and Allocation
//--------------------------------------------------------------
void spectrum_array::init(locate_array tg, locate_array wg,
			  locate_array mg, locate_array pg)
{
  // initialize locate arrays by swapping with the inputs
  time_grid.swap(tg);
  wave_grid.swap(wg);
  mu_grid.swap(mg);
  phi_grid.swap(pg);

  int n_times  = time_grid.size();
  int n_wave   = wave_grid.size();
  int n_mu     = mu_grid.size();
  int n_phi    = phi_grid.size();

  // index parameters
  n_elements  = n_times*n_wave*n_mu*n_phi;
  a3 = n_phi;
  a2 = n_mu*a3;
  a1 = n_wave*a2;

  // allocate memory
  click.resize(n_elements);
  flux.resize(n_elements);

  // clear 
  wipe();
}


//--------------------------------------------------------------
// Functional procedure: Wipe
//--------------------------------------------------------------
void spectrum_array::wipe()
{
  #pragma omp parallel for
  for (int i=0;i<click.size();i++) 
  {
    flux[i]   = 0;
    click[i]  = 0;
  }
}



//--------------------------------------------------------------
// handles the indexing: should be called in this order
//    time, wavelength, mu, phi
//--------------------------------------------------------------
int spectrum_array::index(int t, int l, int m, int p)
{
  return t*a1 + l*a2 + m*a3 + p;
}


//--------------------------------------------------------------
// count a particle
////--------------------------------------------------------------
void spectrum_array::count(double t, double w, double E, double *D)
{
  double mu  = D[2];
  double phi = atan2(D[1],D[0]);

  // if off the LEFT of time/mu/phi grids, just return without counting
  if ((t<time_grid.min) || (mu<mu_grid.min) || (phi<phi_grid.min)) {cout << 1; return;}

  // locate bin number in all dimensions.
  int t_bin = time_grid.locate(t);
  int l_bin = wave_grid.locate(w);
  int m_bin =   mu_grid.locate(mu);
  int p_bin =  phi_grid.locate(phi);

  // if off the RIGHT of time/mu/phi grids, just return without counting
  if((t_bin == time_grid.size()) ||
     (m_bin ==   mu_grid.size()) ||
     (p_bin ==  phi_grid.size())) {cout <<2; return;}

  // if off RIGHT of wavelength grid, store in last bin (LEFT is accounted for by locate)
  if (l_bin == wave_grid.size()) l_bin--;

  // add to counters
  int ind = index(t_bin,l_bin,m_bin,p_bin);

  #pragma omp atomic
  flux[ind]  += E;
  #pragma omp atomic
  click[ind] += 1;
}



//--------------------------------------------------------------
// print out
//--------------------------------------------------------------
void spectrum_array::print()
{
  int nprocs, myID;
  MPI_Comm_size( MPI_COMM_WORLD, &nprocs );
  MPI_Comm_rank( MPI_COMM_WORLD, &myID   );

  if(myID==0){
    FILE *out = fopen(name,"w");

    int n_times  = time_grid.size();
    int n_wave   = wave_grid.size();
    int n_mu     = mu_grid.size();
    int n_phi    = phi_grid.size();

    fprintf(out,"# %d %d %d %d\n",n_times,n_wave,n_mu,n_phi);

    for (int k=0;k<n_mu;k++)
      for (int m=0;m<n_phi;m++)
	for (int i=0;i<n_times;i++)
	  for (int j=0;j<n_wave;j++)
	    {
	      int id = index(i,j,k,m);
	      if (n_times > 1)  fprintf(out,"%12.4e ",time_grid.center(i));;
	      if (n_wave > 1)   fprintf(out,"%12.4e ",wave_grid.center(j));
	      if (n_mu > 1)     fprintf(out,"%12.4f ",mu_grid.center(k));
	      if (n_phi> 1)     fprintf(out,"%12.4f ",phi_grid.center(m));

	      // the delta is infinity if the bin is a catch-all.
	      // Use normalization of 1 to match the hard-coded choice of dt=1 for iterative calculations
	      double wdel = wave_grid.delta(j);
	      double tdel = time_grid.delta(i);
	      double norm = n_mu*n_phi
		* ( wdel < numeric_limits<double>::infinity() ? wdel : 1 )
		* ( tdel < numeric_limits<double>::infinity() ? tdel : 1 );
	      fprintf(out,"%12.5e %10d\n", flux[id]/norm,click[id]);
	    }
    fclose(out);
  }
}


void  spectrum_array::rescale(double r)
{
  #pragma omp parallel for
  for (int i=0;i<flux.size();i++) flux[i] *= r;
}


//--------------------------------------------------------------
// MPI average the spectrum contents
//--------------------------------------------------------------
// only process 0 gets the reduced spectrum to print
void spectrum_array::MPI_average()
{
  int receiving_ID = 0;
  int mpi_procs, myID;

  // average the flux (receive goes out of scope after section)
  {
    vector<double> receive(n_elements,-1);
    MPI_Reduce(&flux.front(), &receive.front(), n_elements, MPI_DOUBLE, MPI_SUM, receiving_ID, MPI_COMM_WORLD);
    flux.swap(receive);
  }

  // average clicks (receive goes out of scope after section)
  {
    vector<int> receive(n_elements,-1);
    MPI_Reduce(&click.front(), &receive.front(), n_elements, MPI_INT, MPI_SUM, receiving_ID, MPI_COMM_WORLD);
    click.swap(receive);
  }

  // only have the receiving ID do the division
  MPI_Comm_size( MPI_COMM_WORLD, &mpi_procs );
  MPI_Comm_rank( MPI_COMM_WORLD, &myID      );
  if(myID == receiving_ID){
    #pragma omp parallel for
    for (int i=0;i<n_elements;i++)
    {
      flux[i]  /= mpi_procs;
      click[i] /= mpi_procs;
    }
  }
}
