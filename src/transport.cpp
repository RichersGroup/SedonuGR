#pragma warning disable 161
#include <omp.h>
#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <iostream>
#include <limits>
#include <cassert>
#include "transport.h"
#include "physical_constants.h"
#include "Lua.h"
#include "grid_1D_sphere.h"
#include "grid_3D_cart.h"
#include "species_general.h"
#include "photons.h"
#include "neutrinos.h"
#include "cdf_array.h"
#include "nulib_interface.h"

namespace pc = physical_constants;

//----------------------------------------------------------------------------
// Initialize the transport module
// Includes setting up the grid, particles,
// and MPI work distribution
//----------------------------------------------------------------------------
void transport::init(Lua* lua)
{ 
  // get mpi rank
  MPI_Comm_size( MPI_COMM_WORLD, &MPI_nprocs );
  MPI_Comm_rank( MPI_COMM_WORLD, &MPI_myID  );
  MPI_real = ( sizeof(real)==4 ? MPI_FLOAT : MPI_DOUBLE );
  verbose = (MPI_myID==0);

  // read simulation parameters
  solve_T       = lua->scalar<int>("solve_T");
  solve_Ye      = lua->scalar<int>("solve_Ye");
  step_size     = lua->scalar<double>("step_size");
  damping       = lua->scalar<double>("damping");
  do_photons    = lua->scalar<int>("do_photons");
  do_neutrinos  = lua->scalar<int>("do_neutrinos");
  if(solve_T || solve_Ye){
    brent_itmax     = lua->scalar<int>("brent_itmax");
    brent_solve_tolerance = lua->scalar<double>("brent_tolerance");
  }

  // figure out what emission models we're using
  n_emit_core  = lua->scalar<int>("n_emit_core");
  n_emit_decay = lua->scalar<int>("n_emit_decay");
  n_emit_therm = lua->scalar<int>("n_emit_therm");
  do_visc      = lua->scalar<int>("do_visc");
  if(do_visc) visc_specific_heat_rate = lua->scalar<double>("visc_specific_heat_rate");

  // Reserve all the memory we might need right now. Speeds up particle additions.
  max_particles = lua->scalar<int>("max_particles");
  particles.reserve(max_particles);

  //=================//
  // SET UP THE GRID //
  //=================//
  // read the grid type
  string grid_type = lua->scalar<string>("grid_type");

  // create a grid of the appropriate type
  if     (grid_type == "grid_1D_sphere") grid = new grid_1D_sphere;
  else if(grid_type == "grid_3D_cart"  ) grid = new grid_3D_cart;
  else{
    if(verbose) std::cout << "ERROR: the requested grid type is not implemented." << std::endl;
    exit(3);}
  
  // initialize the grid (including reading the model file)
  grid->init(lua);

  // calculate integrated quantities to check
  double mass = 0.0;
  double KE   = 0.0;
  #pragma omp parallel for reduction(+:mass,KE)
  for (int i=0;i<grid->z.size();i++){
    double my_mass = grid->z[i].rho * grid->zone_volume(i);
    double* v = grid->z[i].v;
    mass      += my_mass;
    KE        += 0.5 * my_mass * (v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
  }
  if (verbose){
    cout << "# mass = " << mass << " g" <<endl;
    cout << "# KE = " << KE << " erg" << endl;
  }

  //===============//
  // GENERAL SETUP //
  //===============//
  // figure out which zones are in this processors work load
  // a processor will do work in range [start,end)
  my_zone_end.resize(MPI_nprocs);
  for(int proc=0; proc<MPI_nprocs; proc++){
    // how much work does this processor do?
    int my_job = (int)(grid->z.size()/(1.0*MPI_nprocs));
    if(my_job < 1) my_job = 1;

    // where does this processor start and stop its work? (only the end needs to be stored)
    int my_zone_start = proc*my_job;
    my_zone_end[proc] = my_zone_start + my_job;

    // make sure last guy finishes it all
    if(proc == MPI_nprocs-1) my_zone_end[proc] = grid->z.size();

    // make sure nobody goes overboard
    if(my_zone_end[proc] >= grid->z.size()) my_zone_end[proc] = grid->z.size();
  }

  // setup and seed random number generator(s)
  rangen.init();

  //==================//
  // SET UP TRANSPORT //
  //==================//
  // start at time 0
  t_now = 0;

  /**********************/
  /**/ if(do_photons) /**/
  /**********************/
  {
    photons* photons_tmp = new photons;
    photons_tmp->init(lua, this);
    species_list.push_back(photons_tmp);
  }
  /************************/
  /**/ if(do_neutrinos) /**/
  /************************/
  {
    // read the fortran module into memory
    if(verbose) cout << "# Initializing NuLib..." << endl;
    string nulib_table = lua->scalar<string>("nulib_table");
    nulib_init(nulib_table);

    // create a species for each in the nulib table
    int num_nut_species = nulib_get_nspecies();
    for(int i=0; i<num_nut_species; i++){
      neutrinos* neutrinos_tmp = new neutrinos;
      neutrinos_tmp->nulibID = i;
      neutrinos_tmp->num_nut_species = num_nut_species;
      neutrinos_tmp->init(lua, this);
      species_list.push_back(neutrinos_tmp);
    }
  }

  // complain if we're not simulating anything
  if(species_list.size() == 0)
  {
    if(verbose) cout << "ERROR: you must simulate at least one species of particle." << endl;
    exit(7);
  }

  // set global min/max values (make range infinite to catch errors)
  T_min   =  numeric_limits<double>::infinity();
  T_max   = -numeric_limits<double>::infinity();
  Ye_min  =  numeric_limits<double>::infinity();
  Ye_max  = -numeric_limits<double>::infinity();
  rho_min =  numeric_limits<double>::infinity();
  rho_max = -numeric_limits<double>::infinity();
  for(int i=0; i<species_list.size(); i++)
  {
    if(species_list[i]->T_min   < T_min  ) T_min   = species_list[i]->T_min;
    if(species_list[i]->T_max   > T_max  ) T_max   = species_list[i]->T_max;
    if(species_list[i]->Ye_min  < Ye_min ) Ye_min  = species_list[i]->Ye_min;
    if(species_list[i]->Ye_max  > Ye_max ) Ye_max  = species_list[i]->Ye_max;
    if(species_list[i]->rho_min < rho_min) rho_min = species_list[i]->rho_min;
    if(species_list[i]->rho_max > rho_max) rho_max = species_list[i]->rho_max;
  }
  if(T_min >= T_max){
    cout << "ERROR: invalid temperature range." << endl;
    exit(14);
  }
  if(Ye_min >= Ye_max && do_neutrinos){
    cout << "ERROR: invalid Ye range." << endl;
    exit(14);
  }

  // initialize all the zone eas variables
  #pragma omp parallel for collapse(2)
  for(int i=0; i<species_list.size(); i++) 
    for(int j=0; j<grid->z.size(); j++)
      species_list[i]->set_eas(j);

  // scatter initial particles in the simulation area
  // don't initialize them if iterative calculation. They all come from the core.
  // int init_particles = lua->scalar<int>("init_particles");
  // if(!iterate) initialize_particles(init_particles);


  //=================//
  // SET UP THE CORE //
  //=================//
  // the core temperature is used only in setting its emis vector
  // so it's looked at only in species::myInit()
  if(n_emit_core > 0){
    r_core      = lua->scalar<double>("r_core");
    L_core      = lua->scalar<double>("L_core");
    core_species_cdf.resize(species_list.size());
    for(int i=0; i<species_list.size(); i++)
      core_species_cdf.set_value(i, species_list[i]->int_core_emis());
    core_species_cdf.normalize();
  }
  else{
    r_core = 0;
    L_core = 0;
  }
}




//------------------------------------------------------------
// take a transport time step 
//------------------------------------------------------------
void transport::step(const double dt)
{
  assert(dt>0);
  
  // calculate the zone eas variables
  #pragma omp parallel for collapse(2)
  for(int i=0; i<species_list.size(); i++) 
    for(int j=0; j<grid->z.size(); j++)
      species_list[i]->set_eas(j);

  // prepare zone quantities for another round of transport
  double net_visc_heating = 0;
  #pragma omp parallel for reduction(+:net_visc_heating)
  for (int i=0;i<grid->z.size();i++)
  {
    // clear the tallies of the radiation quantities in each zone
    grid->z[i].e_rad    = 0;
    grid->z[i].f_rad[0] = 0;
    grid->z[i].f_rad[1] = 0;
    grid->z[i].f_rad[2] = 0;
    grid->z[i].l_abs    = 0;
    grid->z[i].e_abs    = 0;

    // set the net luminosity (to be recalculated by emission routines)
    L_net = 0;

    // add heat absorbed from viscosity to tally of e_abs
    grid->z[i].e_abs += dt * zone_visc_heat_rate(i);
    net_visc_heating += zone_visc_heat_rate(i);
  }
  if(verbose) cout << "Viscous heating: " << net_visc_heating << " erg/s" << endl;

  // emit new particles
  emit_particles(dt);

  // Propagate the particles
  propagate_particles(dt);

  // properly normalize the radiative quantities
  #pragma omp parallel for
  for (int i=0;i<grid->z.size();i++) 
  {
    double vol = grid->zone_volume(i);
    grid->z[i].e_rad    /= vol*pc::c*dt;
    grid->z[i].e_abs    /= vol*dt;
    grid->z[i].f_rad[0] /= vol*pc::c*dt;
    grid->z[i].f_rad[1] /= vol*pc::c*dt;
    grid->z[i].f_rad[2] /= vol*pc::c*dt;
    grid->z[i].l_abs    /= vol*dt;
  }

  // MPI reduce the radiation quantities so each processor has the information needed to solve its grid values
  if(MPI_nprocs>1) reduce_radiation();

  // solve for T_gas and Ye structure if radiative eq. applied
  if(iterate) solve_eq_zone_values();

  // MPI broadcast the results so all processors have matching fluid properties
  if(MPI_nprocs>1) synchronize_gas();

  // advance time step
  // if (!iterate) t_now += dt;
}




//----------------------------------------------------------------------------
// sum up the number of particles in all species
//----------------------------------------------------------------------------
int transport::total_particles() const{
  return particles.size();
}


//----------------------------------------------------------------------------
// randomly sample the nu-integrated emissivities of all
// species to determine the species of a new particle
// emitted from the core
//----------------------------------------------------------------------------
int transport::sample_core_species() const
{
  // randomly sample the species (precomputed CDF)
  double z = rangen.uniform();
  return core_species_cdf.sample(z);
}



//----------------------------------------------------------------------------
// randomly sample the nu-integrated emissivities of all
// species to determine the species of a new particle
// emitted from a zone
//----------------------------------------------------------------------------
// note: could store a zone_species_cdf structure in transport,
// but this would use more memory. Here, trading CPU cycles for 
// memory. If we are CPU limited, we could change this
int transport::sample_zone_species(const int zone_index) const
{
  cdf_array species_cdf;
  double integrated_emis;
  species_cdf.resize(species_list.size());

  // set values and normalize
  for(int i=0; i<species_list.size(); i++)
  {
    integrated_emis = species_list[i]->int_zone_emis(zone_index);
    species_cdf.set_value(i,integrated_emis);
  }
  species_cdf.normalize();

  // randomly sample the species
  double z = rangen.uniform();
  return species_cdf.sample(z);
}


//------------------------------------------------------------
// Combine the radiation tallies in all zones
// from all processors using MPI allreduce
//------------------------------------------------------------
void transport::reduce_radiation()
{
  vector<real> send, receive;
  int my_begin, my_end, size;

  //-- EACH PROCESSOR GETS THE REDUCTION INFORMATION IT NEEDS
  for(int proc=0; proc<MPI_nprocs; proc++){

    // set the begin and end indices so a process covers range [begin,end)
    my_begin = ( proc==0 ? 0 : my_zone_end[proc-1] );
    my_end = my_zone_end[proc];

    // set the computation size and create the send/receive vectors
    size = my_end - my_begin;
    send.resize(size);
    receive.resize(size);

    // reduce e_abs
    if(solve_T){
      for(int i=my_begin; i<my_end; i++) send[i-my_begin] = grid->z[i].e_abs;
      MPI_Reduce(&send.front(), &receive.front(), size, MPI_real, MPI_SUM, proc, MPI_COMM_WORLD);
      for(int i=my_begin; i<my_end; i++) grid->z[i].e_abs = receive[i-my_begin] / (real)MPI_nprocs;
    }

    // reduce l_abs
    if(solve_Ye){
      for(int i=my_begin; i<my_end; i++) send[i-my_begin] = grid->z[i].l_abs;
      MPI_Reduce(&send.front(), &receive.front(), size, MPI_real, MPI_SUM, proc, MPI_COMM_WORLD);
      for(int i=my_begin; i<my_end; i++) grid->z[i].l_abs = receive[i-my_begin] / (real)MPI_nprocs;
    }

    // TODO - need to put in other quantities...
  }
}

void transport::synchronize_gas()
{
  vector<real> buffer;
  int my_begin, my_end, size;

  //-- EACH PROCESSOR SENDS THE GRID INFORMATION IT SOLVED
  for(int proc=0; proc<MPI_nprocs; proc++){

    // set the begin and end indices so a process covers range [begin,end)
    my_begin = ( proc==0 ? 0 : my_zone_end[proc-1] );
    my_end = my_zone_end[proc];

    // set the computation size and create the send/receive vectors
    size = my_end - my_begin;
    buffer.resize(size);

    // broadcast T_gas
    if(solve_T){
      if(proc==MPI_myID) for(int i=my_begin; i<my_end; i++) buffer[i-my_begin] = grid->z[i].T_gas;
      MPI_Bcast(&buffer.front(), size, MPI_real, proc, MPI_COMM_WORLD);
      if(proc!=MPI_myID) for(int i=my_begin; i<my_end; i++) grid->z[i].T_gas = buffer[i-my_begin];
    }

    // broadcast Ye
    if(solve_Ye){
      if(proc==MPI_myID) for(int i=my_begin; i<my_end; i++) buffer[i-my_begin] = grid->z[i].Ye;
      MPI_Bcast(&buffer.front(), size, MPI_real, proc, MPI_COMM_WORLD);
      if(proc!=MPI_myID) for(int i=my_begin; i<my_end; i++) grid->z[i].Ye = buffer[i-my_begin];
    }
  }
}


// rate at which viscosity energizes the fluid (erg/s)
double transport::zone_visc_heat_rate(const int zone_index) const{
  return visc_specific_heat_rate * grid->z[zone_index].rho * grid->zone_volume(zone_index);
}


// Write the spectra to disk (using data since the last write)
void transport::write_spectra(const int it)
{
  for(int i=0; i<species_list.size(); i++){
    if(MPI_nprocs>1) species_list[i]->spectrum.MPI_average();
    char sname[100];
    sprintf(sname,"species%d_I%d.spec",i,it);
    species_list[i]->spectrum.set_name(sname);
    species_list[i]->spectrum.print();
    species_list[i]->spectrum.wipe();
  }
}
