#include <mpi.h>
#include <gsl/gsl_rng.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <iomanip>
#include <math.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include "Lua.h"
#include "grid_1D_sphere.h"
#include "transport.h"
#include "physical_constants.h"

namespace pc = physical_constants;

//--------------------------------------------------------
// The main code
// The user writes this for their own needs
//--------------------------------------------------------
int main(int argc, char **argv)
{
  int read_model_file(string model_file, grid_general**);
  void write_zones(int, grid_general *this_grid);

  // set the global start timer
  time_t start_tp,end_tp;
  time(&start_tp);

  // initialize MPI parallelism
  int my_rank,n_procs;
  MPI_Init( &argc, &argv );
  MPI_Comm_rank( MPI_COMM_WORLD, &my_rank );
  MPI_Comm_size( MPI_COMM_WORLD, &n_procs);
  
  // verbocity
  const int verbose = (my_rank == 0);
  if (verbose) cout << "# MPI cores used = " << n_procs << endl;
  
  // start timer 
  double proc_time_start = MPI_Wtime();

  // open up the lua parameter file
  Lua lua;
  std::string    script_file = "param.lua";
  if( argc > 1 ) script_file = std::string( argv[ 1 ] );
  lua.init( script_file );

  // read initial setup from modelfile
  grid_general *grid;
  std::string model_file = lua.scalar<std::string>("model_file"); 
  read_model_file(model_file,&grid);

  // calculate integrated quantities to check
  double mass = 0.0;
  double KE   = 0.0;
  for (int i=0;i<grid->n_zones;i++) 
  {
    mass += grid->z[i].rho*grid->zone_volume(i);
  }
  if (verbose) cout << "# mass = " << mass << endl;

  // read in time stepping parameters
  int    n_times     = lua.scalar<int>("n_times");
  int    iterate     = lua.scalar<int>("iterate");

  double t_stop      = lua.scalar<double>("t_stop");
  double tstep_max   = lua.scalar<double>("tstep_max");
  double tstep_min   = lua.scalar<double>("tstep_min");
  double tstep_start = lua.scalar<double>("tstep_start");
  double tstep_del   = lua.scalar<double>("tstep_del");
  double write_out   = lua.scalar<double>("write_out");

  // set up transport module
  int use_transport  = lua.scalar<int>("use_transport");
  transport transport;
  transport.init(script_file,grid);
  
  // loop over time steps;
  double t = 0;
  double t_step = tstep_start;
  int   iw = 0; // number of times output has been written
 
  // check for a time static, iterative calculation
  if (iterate) n_times = iterate;

  // loop over times
  for (int it=0;it<n_times;it++)
  {
    // get this time step (ignored if iterative calc)
    if (t_step < tstep_min) t_step = tstep_min;
    if (t_step > tstep_max) t_step = tstep_max;
    if ((tstep_del > 0)&&(t > 0)) 
      if (t_step > t*tstep_del) t_step = t*tstep_del;

    // printout time step
    if (verbose)
      printf("%8d %12.4e %12.4e %5d\n",it,t,t_step, transport.num_particles());

    // writeout zone state when appropriate 
    if ((verbose)&&((t >= write_out*iw)||(iterate)))
    {
      printf("# writing zone %d at time %e\n",iw, t);
      write_zones(iw,grid);
      iw++;
    }

    // do transport step
    transport.step(t_step);

    // print out spectrum in iterative calc
    if (iterate)
    {
      char sname[100];
      sprintf(sname,"optical_I%d.spec",it+1);
      transport.spectrum.set_name(sname);
      transport.spectrum.MPI_average();
      transport.spectrum.print();
      transport.spectrum.wipe();
    }

    // advance time
    if (!iterate) t = t + t_step;
    if (t > t_stop) break;
  }

  // print out final spectrum
  if (!iterate)  {
    transport.spectrum.MPI_average();
    transport.spectrum.print(); }

  //---------------------------------------------------------------------
  // CALCULATION DONE; WRITE OUT AND FINISH
  //---------------------------------------------------------------------

  // calculate the elapsed time 
  double proc_time_end = MPI_Wtime();
  double time_wasted = proc_time_end - proc_time_start;
  
  if (verbose)
    printf("#\n# CALCULATION took %.3e seconds or %.3f mins or %.3f hours\n",
	   time_wasted,time_wasted/60.0,time_wasted/60.0/60.0);

  // finish up mpi
  MPI_Finalize();


  return 0;
}

//---------------------------------------------------------------------
// user writes their own function that fills in the initial
// state of the grid.  In this case, we read in from a file
//---------------------------------------------------------------------
int read_model_file(string model_file, grid_general **grid)
{
  std::ifstream infile;   

  infile.open(model_file.c_str());

  // check if open worked
  if(infile.fail()) 
  { 
    cout << "error, can't read model file" << model_file << endl; 
    return 1; 
  } 


  // geometry of model
  string grid_type;
  infile >> grid_type;

  // ********************************************
  // read in and set up 1D model
  // ********************************************
  if (grid_type == "1D_sphere") 
  {
    // allocate model
    grid_1D_sphere *this_grid;
    this_grid = new grid_1D_sphere;

    // number of zones
    int nx;
    infile >> nx;

    // read radii of zones
    double r_in;
    std::vector<double> r_out, dens, ni, T;
    double x;
    infile >> r_in;
    
    double texp;
    infile >> texp;

    for (int i=0;i<nx; i++)
    { 
      infile >> x;
      r_out.push_back(x);
      infile >> x;
      dens.push_back(x);
      infile >> x;
      T.push_back(x);
      infile >> x;
      ni.push_back(x);
    } 

    // initialize the grid
    this_grid->init(r_in,r_out);

    // fill up zones 
    for (int i=0;i<nx;i++) 
    {
      this_grid->z[i].rho    = dens[i];
      this_grid->z[i].ni56   = ni[i]; 
      this_grid->z[i].v[0]   = r_out[i]/texp;
      this_grid->z[i].T_gas  = T[i];
      this_grid->z[i].e_rad  = pc::a*pow(T[i],4);
    }
  
    // set the main grid to this data
    *grid = this_grid;
  }
  
  infile.close(); 

  return 0;
}


void write_zones(int iw, grid_general *this_grid)
{
  char zonefile[1000];
  char base[1000];

  if (iw < 10) sprintf(base,"_0000%d",iw);
  else if (iw < 100) sprintf(base,"_000%d",iw);
  else if (iw < 1000) sprintf(base,"_00%d",iw);
  else if (iw < 10000) sprintf(base,"_0%d",iw);
  else sprintf(base,"_%d",iw);
  sprintf(zonefile,"ray%s",base);

  ofstream outf;
  outf.open(zonefile);
  //  outf << setw(12);
  outf << setprecision(4);
  outf << scientific;

  for (int i=0;i<this_grid->n_zones;i++)
  {
    double r[3];
    this_grid->coordinates(i,r); 
    outf << r[0] << " ";
    zone z = this_grid->z[i];

    double T_rad = pow(z.e_rad/pc::a,0.25);
    outf << T_rad << " ";
    outf << z.T_gas << " ";
    
    outf << endl;
  }

}
