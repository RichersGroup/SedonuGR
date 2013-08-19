#include "transport.h"
#include <limits>
#include <mpi.h>
#include <gsl/gsl_rng.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <iostream>
#include "physical_constants.h"
#include "Lua.h"
#include "grid_1D_sphere.h"
#include "grid_3D_cart.h"
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
  int my_rank, n_procs;
  MPI_Comm_size( MPI_COMM_WORLD, &n_procs );
  MPI_Comm_rank( MPI_COMM_WORLD, &my_rank );
  verbose = (my_rank==0);

  //=================//
  // SET UP THE GRID //
  //=================//
  // read the grid type
  string grid_type = lua->scalar<string>("grid_type");

  // create a grid of the appropriate type
  if     (grid_type == "grid_1D_sphere") grid = new grid_1D_sphere;
  else if(grid_type == "grid_3D_cart"  ) grid = new grid_3D_cart;
  else{
    if(verbose) std::cout << "Error: the requested grid type is not implemented." << std::endl;
    exit(3);}
  
  // initialize the grid (including reading the model file)
  grid->init(lua);

  // calculate integrated quantities to check
  double mass = 0.0;
  double KE   = 0.0;
  for (int i=0;i<grid->z.size();i++)
  {
    mass += grid->z[i].rho * grid->zone_volume(i);
  }
  if (verbose) cout << "# mass = " << mass << endl;

  //===============//
  // GENERAL SETUP //
  //===============//
  // figure out which zones are in this processors work load
  int my_job = (int)(grid->z.size()/(1.0*n_procs));
  if (my_job < 1) my_job = 1;
  my_zone_start = my_rank*my_job;
  my_zone_end   = my_zone_start + my_job;
  // make sure last guy finishes it all
  if (my_rank == n_procs-1) my_zone_end = grid->z.size();
  // get rid of unneeded processers
  if (my_zone_start >= grid->z.size()) my_zone_start = my_zone_end;

  // setup and seed random number generator
  const gsl_rng_type * TypeR;
  gsl_rng_env_setup();
  gsl_rng_default_seed = (unsigned int)time(NULL) + my_rank;
  TypeR = gsl_rng_default;
  rangen = gsl_rng_alloc (TypeR);


  //==================//
  // SET UP TRANSPORT //
  //==================//
  // start at time 0
  t_now = 0;

  // read relevant parameters
  radiative_eq  = lua->scalar<int>("radiative_eq");
  iterate       = lua->scalar<int>("iterate");
  step_size     = lua->scalar<double>("step_size");

  // determine which species to simulate
  int do_photons   = lua->scalar<int>("do_photons");
  int do_neutrinos = lua->scalar<int>("do_neutrinos");

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
    cout << "Initializing NuLib..." << endl;
    string nulib_table = lua->scalar<string>("nulib_table");
    nulib_init(nulib_table);
    neutrinos* neutrinos_tmp;

    // create a species for each in the nulib table
    int num_nut_species = nulib_get_nspecies();
    for(int i=0; i<num_nut_species; i++){
      neutrinos_tmp = new neutrinos;
      neutrinos_tmp->nulibID = i;
      neutrinos_tmp->init(lua, this);
      species_list.push_back(neutrinos_tmp);
    }
  }

  // complain if we're not simulating anything
  if(species_list.size() == 0)
  {
    if(verbose) cout << "Error: you must simulate at least one species of particle." << endl;
    exit(7);
  }

  // initialize all the zone eas variables
  for(int i=0; i<species_list.size(); i++) 
    for(int j=0; j<species_list[i]->size(); j++)
      species_list[i]->set_eas(j);

  // scatter initial particles in the simulation area
  int init_particles = lua->scalar<int>("init_particles");
  initialize_particles(init_particles);

  //=================//
  // SET UP THE CORE //
  //=================//
  // the core temperature is used only in setting its emis vector
  // so it's looked at only in species::myInit()
  n_inject = lua->scalar<int>("n_inject");
  r_core   = lua->scalar<double>("r_core");
  L_core   = lua->scalar<double>("L_core");
  core_species_cdf.resize(species_list.size());
  for(int i=0; i<species_list.size(); i++){
    core_species_cdf.set_value(i, species_list[i]->int_core_emis());}
  core_species_cdf.normalize();
}




//------------------------------------------------------------
// take a transport time step 
//------------------------------------------------------------
void transport::step(double dt)
{
  int i;

  // nominal time for iterative calc is 1
  if (iterate) dt = 1;
  
  // calculate the zone eas variables
  for(int i=0; i<species_list.size(); i++) 
    for(int j=0; j<species_list[i]->size(); j++)
      species_list[i]->set_eas(j);

  // emit new particles
  emit_particles(dt);

  // clear the tallies of the radiation quantities in each zone
  for (int i=0;i<grid->n_zones;i++) 
  {
    grid->z[i].e_rad  = 0;
    grid->z[i].e_rad  = 0;
    grid->z[i].e_abs  = 0;
    grid->z[i].fx_rad = 0;
    grid->z[i].fy_rad = 0;
    grid->z[i].fz_rad = 0;
  }

  // Propagate the particles
  for(int i=0; i<species_list.size(); i++) species_list[i]->propagate_particles(dt);

  // properly normalize the radiative quantities
  for (int i=0;i<grid->n_zones;i++) 
  {
    double vol = grid->zone_volume(i);
    grid->z[i].e_rad   /= vol*pc::c*dt;
    grid->z[i].e_abs   /= vol*dt; 
    grid->z[i].fx_rad  /= vol*pc::c*dt; 
    grid->z[i].fy_rad  /= vol*pc::c*dt;
    grid->z[i].fz_rad  /= vol*pc::c*dt;
  }

  // MPI reduce the tallies and put in place
  grid->reduce_radiation();

  // solve for T_gas structure if radiative eq. applied
  if (radiative_eq) solve_eq_temperature();
   
  // advance time step
  if (!iterate) t_now += dt;
}



//-------------------------------------------------------------
//  Solve for the temperature assuming radiative equilibrium
//-------------------------------------------------------------
void transport::solve_eq_temperature()
{
  // wipe temperature structure
  for (int i=0;i<grid->n_zones;i++) grid->z[i].T_gas = 0;

  // solve radiative equilibrium temperature
  for (int i=this->my_zone_start;i<this->my_zone_end;i++)
    grid->z[i].T_gas = temp_brent_method(i);

  // mpi reduce the results
  grid->reduce_T_gas();
}


//----------------------------------------------------------------------------
// This is the function that expresses radiative equillibrium
// in a cell (i.e. E_absorbed = E_emitted).  It is used in
// The Brent solver below to determine the temperature such
// that RadEq holds
//----------------------------------------------------------------------------
double transport::rad_eq_function(int zone_index,double T)
{
  // total energy absorbed in zone
  double E_absorbed = grid->z[zone_index].e_abs;
  // total energy emitted (to be calculated)
  double E_emitted = 0.;

  // set the zone temperature
  grid->z[zone_index].T_gas = T;
  
  // include the emission from all species
  for(int i=0; i<species_list.size(); i++)
  {
    // reset the eas variables in this zone
    // OPTIMIZE - only set the emissivity variable
    species_list[i]->set_eas(zone_index);

    // integrate emisison over frequency (angle
    // integration gives the 4*PI) to get total
    // radiation energy emitted. Opacities are
    // held constant for this (assumed not to change
    // much from the last time step).
    E_emitted += 4.0*pc::pi * species_list[i]->int_zone_emis(zone_index);
  }

  // radiative equillibrium condition: "emission equals absorbtion"
  // return to Brent function to iterate this to zero
  return (E_emitted - E_absorbed);
}


//-----------------------------------------------------------
// Brents method (from Numerical Recipes) to solve 
// non-linear equation for T in rad equillibrium
//-----------------------------------------------------------
// definitions used for temperature solver
#define TEMP_SOLVE_TOLERANCE 1.e-2
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
double transport::temp_brent_method(int cell)
{  
  int ITMAX = 100;
  double EPS = 3.0e-8;
  int iter;

  // Initial guesses
  double a=TEMP_RANGE_MIN;
  double b=TEMP_RANGE_MAX;
  double c=b;
  double d,e,min1,min2;
  double fa=rad_eq_function(cell,a);
  double fb=rad_eq_function(cell,b);
  double fc,p,q,r,s,tol1,xm;
  
  //if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0))
  //  printf("Root must be bracketed in zbrent");
  fc=fb;
  for (iter=1;iter<=ITMAX;iter++) {
    if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
      c=a;
      fc=fa;
      e=d=b-a;
    }
    if (fabs(fc) < fabs(fb)) {
      a=b;
      b=c;
      c=a;
      fa=fb;
      fb=fc;
      fc=fa;
    }
    tol1=2.0*EPS*fabs(b)+0.5*TEMP_SOLVE_TOLERANCE;
    xm=0.5*(c-b);
    if (fabs(xm) <= tol1 || fb == 0.0) return b;
    if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
      s=fb/fa;
      if (a == c) {
         p=2.0*xm*s;
         q=1.0-s;
      } else {
	q=fa/fc;
	r=fb/fc;
	p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
	q=(q-1.0)*(r-1.0)*(s-1.0);
      }
      if (p > 0.0) q = -q;
      p=fabs(p);
      min1=3.0*xm*q-fabs(tol1*q);
      min2=fabs(e*q);
      if (2.0*p < (min1 < min2 ? min1 : min2)) {
	e=d;
	d=p/q;
      } else {
	d=xm;
	e=d;
      }
    } else {
      d=xm;
      e=d;
    }
    a=b;
    fa=fb;
    if (fabs(d) > tol1)
      b += d;
    else
      b += SIGN(tol1,xm);
    fb=rad_eq_function(cell,b);
  }
  printf("Maximum number of iterations exceeded in zbrent");
  return 0.0;
}

//----------------------------------------------------------------------------
// sum up the number of particles in all species
//----------------------------------------------------------------------------
int transport::total_particles()
{
  int total = 0;
  for(int i=0; i<species_list.size(); i++) total += species_list[i]->size();
  return total;
}


//----------------------------------------------------------------------------
// randomly sample the nu-integrated emissivities of all
// species to determine the species of a new particle
// emitted from the core
//----------------------------------------------------------------------------
int transport::sample_core_species()
{
  // randomly sample the species (precomputed)
  double z = gsl_rng_uniform(rangen);
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
int transport::sample_zone_species(int zone_index)
{
  cdf_array species_cdf;
  double integrated_emis;

  // set values and normalize
  for(int i=0; i<species_list.size(); i++)
  {
    integrated_emis = species_list[i]->int_zone_emis(zone_index);
    species_cdf.set_value(i,integrated_emis);
  }
  species_cdf.normalize();

  // randomly sample the species
  double z = gsl_rng_uniform(rangen);
  return species_cdf.sample(z);
}
