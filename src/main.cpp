#include <mpi.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <iomanip>
#include <math.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include "Lua.h"
#include "transport.h"
#include "species_general.h"
#include "physical_constants.h"


namespace pc = physical_constants;
using namespace std;

//--------------------------------------------------------
// The main code
// The user writes this for their own needs
//--------------------------------------------------------
int main(int argc, char **argv)
{
  //============//
  // INITIALIZE //
  //============//
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

  // set up the transport module (includes the grid)
  cout << "initializing" << endl;
  transport sim;
  sim.init(&lua);

  // read in time stepping parameters
  int    n_times     = lua.scalar<int>("n_times");
  int    iterate     = lua.scalar<int>("iterate");
  double t_stop      = lua.scalar<double>("t_stop");
  double tstep_max   = lua.scalar<double>("tstep_max");
  double tstep_min   = lua.scalar<double>("tstep_min");
  double tstep_start = lua.scalar<double>("tstep_start");
  double tstep_del   = lua.scalar<double>("tstep_del");
  double write_out   = lua.scalar<double>("write_out");
  lua.close();

  // loop over time steps;
  double t = 0;
  double t_step = tstep_start;
  int   iw = 0; // number of times output has been written

  // check for a time static, iterative calculation
  if (iterate) n_times = iterate;


  //===========//
  // TIME LOOP //
  //===========//
  for (int it=0;it<n_times;it++)
  {
    // get this time step (ignored if iterative calc)
    if (t_step < tstep_min) t_step = tstep_min;
    if (t_step > tstep_max) t_step = tstep_max;
    if ( (tstep_del>0) && (t>0) && (t_step>t*tstep_del) ) t_step = t*tstep_del;

    // printout time step
    if (verbose)
      printf("%8d %12.4e %12.4e %5d\n",it,t,t_step, sim.total_particles());

    // writeout zone state when appropriate 
    if ( verbose && ( (t>=write_out*iw) || (iterate) ) )
    {
      printf("# writing zone %d at time %e\n",iw, t);
      sim.grid->write_zones(iw);
      iw++;
    }

    // do transport step
    sim.step(t_step);

    // print out spectrum in iterative calc
    if (iterate) for(int i=0; i<sim.species_list.size(); i++)
    {
      char sname[100];
      sprintf(sname,"optical_I%d.spec",it+1);
      sim.species_list[i]->spectrum.set_name(sname);
      sim.species_list[i]->spectrum.MPI_average();
      sim.species_list[i]->spectrum.print();
      sim.species_list[i]->spectrum.wipe();
    }

    // advance time
    if (!iterate) t = t + t_step;
    if (t > t_stop) break;
  }

  //===============================//
  // PRINT FINAL SPECTRUM AND EXIT //
  //===============================//
  // print final spectrum
  if (!iterate) for(int i=0; i<sim.species_list.size(); i++){
    sim.species_list[i]->spectrum.MPI_average();
    sim.species_list[i]->spectrum.print(); 
  }

  // calculate the elapsed time 
  double proc_time_end = MPI_Wtime();
  double time_wasted = proc_time_end - proc_time_start;
  if (verbose)
    printf("#\n# CALCULATION took %.3e seconds or %.3f mins or %.3f hours\n",
	   time_wasted,time_wasted/60.0,time_wasted/60.0/60.0);

  // exit the program
  delete sim.grid;
  for(int i=0; i<sim.species_list.size(); i++) delete sim.species_list[i];
  MPI_Finalize();
  return 0;
}
