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
  do_photons   = lua->scalar<int>("do_photons");
  do_neutrinos = lua->scalar<int>("do_neutrinos");

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
      neutrinos_tmp->num_nut_species = num_nut_species;
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

  // set global min/max values
  T_min  = species_list[0]->T_min;
  T_max  = species_list[0]->T_max;
  Ye_min = species_list[0]->Ye_min;
  Ye_max = species_list[0]->Ye_max;
  for(int i=0; i<species_list.size(); i++)
  {
    if(species_list[i]->T_min < T_min) T_min = species_list[i]->T_min;
    if(species_list[i]->T_max > T_max) T_max = species_list[i]->T_max;
    if(species_list[i]->Ye_min < Ye_min) Ye_min = species_list[i]->Ye_min;
    if(species_list[i]->Ye_max > Ye_max) Ye_max = species_list[i]->Ye_max;
  }
  if(T_min >= T_max){
    cout << "ERROR: invalid temperature range." << endl;
    exit(14);
  }
  if(Ye_min >= Ye_max && do_neutrinos){
    cout << "ERROR: invalid Ye range." << endl;
    exit(14);
  }
  cout << "global T_min is " << T_min << endl;
  cout << "global T_max is " << T_max << endl;

  // initialize all the zone eas variables
  for(int i=0; i<species_list.size(); i++) 
    for(int j=0; j<grid->z.size(); j++)
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
    for(int j=0; j<grid->z.size(); j++)
      species_list[i]->set_eas(j);

  // emit new particles
  emit_particles(dt);

  // clear the tallies of the radiation quantities in each zone
  for (int i=0;i<grid->z.size();i++) 
  {
    grid->z[i].e_rad  = 0;
    grid->z[i].e_rad  = 0;
    grid->z[i].e_abs  = 0;
    grid->z[i].fx_rad = 0;
    grid->z[i].fy_rad = 0;
    grid->z[i].fz_rad = 0;
    grid->z[i].l_abs  = 0;
  }

  // Propagate the particles
  for(int i=0; i<species_list.size(); i++) species_list[i]->propagate_particles(dt);

  // properly normalize the radiative quantities
  for (int i=0;i<grid->z.size();i++) 
  {
    double vol = grid->zone_volume(i);
    grid->z[i].e_rad   /= vol*pc::c*dt;
    grid->z[i].e_abs   /= vol*dt; 
    grid->z[i].fx_rad  /= vol*pc::c*dt; 
    grid->z[i].fy_rad  /= vol*pc::c*dt;
    grid->z[i].fz_rad  /= vol*pc::c*dt;
    grid->z[i].l_abs   /= vol*dt;
  }

  // MPI reduce the tallies and put in place
  grid->reduce_radiation();

  // solve for T_gas structure if radiative eq. applied
  if (radiative_eq) solve_eq_zone_values();
   
  // apply changes to composition
  //update_composition();

  // advance time step
  if (!iterate) t_now += dt;
}



//----------------------------------------------------------------------------
// This is the function that expresses radiative equillibrium
// in a cell (i.e. E_absorbed = E_emitted).  It is used in
// The Brent solver below to determine the temperature such
// that RadEq holds
//----------------------------------------------------------------------------
double temp_eq_function(int zone_index, double T, transport* sim)
{
  // total energy absorbed in zone
  double E_absorbed = sim->grid->z[zone_index].e_abs;
  // total energy emitted (to be calculated)
  double E_emitted = 0.;

  // set the zone temperature
  sim->grid->z[zone_index].T_gas = T;
  
  // include the emission from all species
  for(int i=0; i<sim->species_list.size(); i++)
  {
    // reset the eas variables in this zone
    // OPTIMIZE - only set the emissivity variable
    sim->species_list[i]->set_eas(zone_index);

    // integrate emisison over frequency (angle
    // integration gives the 4*PI) to get total
    // radiation energy emitted. Opacities are
    // held constant for this (assumed not to change
    // much from the last time step).
    E_emitted += 4.0*pc::pi * sim->species_list[i]->int_zone_emis(zone_index);
  }
  
  // radiative equillibrium condition: "emission equals absorbtion"
  // return to Brent function to iterate this to zero
  return (E_emitted - E_absorbed);
}


//----------------------------------------------------------------------------
// This is the function that expresses radiative equillibrium
// in a cell (i.e. E_absorbed = E_emitted).  It is used in
// The Brent solver below to determine the temperature such
// that RadEq holds
//----------------------------------------------------------------------------
double Ye_eq_function(int zone_index, double Ye, transport* sim)
{
  // total energy absorbed in zone
  double l_absorbed = sim->grid->z[zone_index].l_abs;
  // total energy emitted (to be calculated)
  double l_emitted = 0.;

  // set the zone temperature
  sim->grid->z[zone_index].Ye = Ye;
  
  // include the emission from all species
  for(int i=0; i<sim->species_list.size(); i++)
  {
    // reset the eas variables in this zone
    // OPTIMIZE - only set the emissivity variable
    sim->species_list[i]->set_eas(zone_index);

    // integrate emisison over frequency (angle
    // integration gives the 4*PI) to get total
    // radiation energy emitted. Opacities are
    // held constant for this (assumed not to change
    // much from the last time step).
    l_emitted += 4.0*pc::pi * sim->species_list[i]->int_zone_lepton_emis(zone_index);
  }

  // radiative equillibrium condition: "emission equals absorbtion"
  // return to Brent function to iterate this to zero
  return (l_emitted - l_absorbed);
}


//-------------------------------------------------------------
//  Solve for the temperature assuming radiative equilibrium
//-------------------------------------------------------------
#define BRENT_SOLVE_TOLERANCE 1.e-2
#define BRENT_ITMAX 100
void transport::solve_eq_zone_values()
{
  double T_last,Ye_last;
  double T_error,Ye_error;
  int iter;
  // solve radiative equilibrium temperature and Ye
  // TODO - loop through this until we reach the brent tolerance for both
  for (int i=my_zone_start; i<my_zone_end; i++)
  {
    iter = 0;
    T_error  = 10*BRENT_SOLVE_TOLERANCE;
    Ye_error = 10*BRENT_SOLVE_TOLERANCE;

    // loop through solving the temperature and Ye until both are within error.
    while(iter<=BRENT_ITMAX && (T_error>BRENT_SOLVE_TOLERANCE || Ye_error>BRENT_SOLVE_TOLERANCE))
    {
      T_last  = grid->z[i].T_gas;
      Ye_last = grid->z[i].Ye;
      grid->z[i].T_gas = brent_method(i, temp_eq_function, T_min,  T_max);
      grid->z[i].Ye    = brent_method(i, Ye_eq_function,   Ye_min, Ye_max);
      T_error  = fabs( (grid->z[i].T_gas - T_last ) / (T_last ) );
      Ye_error = fabs( (grid->z[i].Ye    - Ye_last) / (Ye_last) );
      iter++;
    }

    // warn if it didn't converge
    if(iter == BRENT_ITMAX) cout << "WARNING: outer Brent solver hit maximum iterations." << endl;
  }
  
  // mpi reduce the results
  grid->reduce_gas();
}




//-----------------------------------------------------------
// Brents method (from Numerical Recipes) to solve 
// non-linear equation for T in rad equillibrium
//-----------------------------------------------------------
// definitions used for temperature solver
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
double transport::brent_method(int zone_index, double (*eq_function)(int,double,transport*), double min, double max)
{
  double small = 3.0e-8;
  int iter;

  // Initial guesses
  double a=min;
  double b=max;
  double c=b;
  double d,e,min1,min2;
  double fa=(*eq_function)(zone_index,a,this);
  double fb=(*eq_function)(zone_index,b,this);
  double fc,p,q,r,s,tol1,xm;
  
  //if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0))
  //  printf("Root must be bracketed in zbrent");
  fc=fb;
  for (iter=1;iter<=BRENT_ITMAX;iter++) {
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
    tol1=2.0*small*fabs(b)+0.5*BRENT_SOLVE_TOLERANCE;
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
    fb=(*eq_function)(zone_index,b,this);
  }
  printf("Maximum number of iterations exceeded in zbrent\n");
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


//--------------------//
// update_composition //
//--------------------//
void transport::update_composition()
{
  double lepton_density;
  for(int i=0; i<grid->z.size(); i++)
  {
    // first find what the lepton density was throughout the timestep
    // ACCURACY - note this assumes that m_n = m_p to avoid subtractive cancellation issues.
    // reduces accuracy of Ye to about .15%
    lepton_density = grid->z[i].rho/pc::m_p * grid->z[i].Ye ;

    // then add l_abs (which was previously divided by the zone volume
    lepton_density += grid->z[i].l_abs;

    // convert back into electron fraction
    grid->z[i].Ye = lepton_density*pc::m_p / grid->z[i].rho;
  }
}
