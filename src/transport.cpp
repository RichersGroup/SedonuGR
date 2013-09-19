#include <omp.h>
#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <iostream>
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
  int my_rank, n_procs;
  MPI_Comm_size( MPI_COMM_WORLD, &n_procs );
  MPI_Comm_rank( MPI_COMM_WORLD, &my_rank );
  verbose = (my_rank==0);

  // read simulation parameters
  radiative_eq = lua->scalar<int>("radiative_eq");
  iterate      = lua->scalar<int>("iterate");
  step_size    = lua->scalar<double>("step_size");
  damping      = lua->scalar<double>("damping");
  solve_T      = lua->scalar<int>("solve_T");
  solve_Ye     = lua->scalar<int>("solve_Ye");
  do_photons   = lua->scalar<int>("do_photons");
  do_neutrinos = lua->scalar<int>("do_neutrinos");

  // figure out what zone emission models we're using
  int n_emit_heat  = lua->scalar<int>("n_emit_heat");
  int n_emit_visc  = lua->scalar<int>("n_emit_visc");
  bool do_heat     = n_emit_heat  > 0;
  bool do_visc     = n_emit_visc  > 0;
  do_therm         = ( radiative_eq ? do_visc     : do_heat     );
  n_emit_therm     = ( radiative_eq ? n_emit_visc : n_emit_heat );

  n_emit_decay     = lua->scalar<int>("n_emit_decay");
  do_decay         = n_emit_decay > 0;

  n_emit_core      = lua->scalar<int>("n_emit_core");
  do_core          = n_emit_core>0;

  // complain if the parameters don't make sense together
  if(n_emit_visc>0) visc_specific_heat_rate = lua->scalar<int>("visc_specific_heat_rate");
  if(do_heat && do_visc){
    cout << "ERROR: n_emit_heat and n_emit_visc cannot both be greater than 0." << endl;
    exit(1);
  }
  if(do_heat && radiative_eq){
    cout << "ERROR: n_emit_heat requires that radiative_eq==0" << endl;
    exit(1);
  }
  if(do_visc && !radiative_eq){
    cout << "ERROR: n_emit_visc requires that radiative_eq==1" << endl;
    exit(1);
  }

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
  double therm_lum=0, decay_lum=0;
  #pragma omp parallel for reduction(+:mass,KE,therm_lum,decay_lum)
  for (int i=0;i<grid->z.size();i++){
    double my_mass = grid->z[i].rho * grid->zone_volume(i);
    double* v = grid->z[i].v;
    mass      += my_mass;
    KE        += 0.5 * my_mass * (v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    if(do_therm) therm_lum += ( radiative_eq ? zone_visc_heat_rate(i) : zone_heat_lum(i) );
    if(do_decay) decay_lum += zone_decay_lum(i);
  }
  if (verbose){
    cout << "# mass = " << mass << " g" <<endl;
    cout << "# KE = " << KE << " erg" << endl;
    if(do_heat)  cout << "# L_heat = "  << (radiative_eq ? 0 : therm_lum)  << " erg/s" << endl;
    if(do_visc)  cout << "# L_visc = "  << (radiative_eq ? therm_lum : 0)  << " erg/s" << endl;
    if(do_decay) cout << "# L_decay = " << decay_lum                       << " erg/s" << endl;
  }

  //===============//
  // GENERAL SETUP //
  //===============//
  // figure out which zones are in this processors work load
  int my_job = (int)(grid->z.size()/(1.0*n_procs));
  if (my_job < 1) my_job = 1;
  my_zone_start = my_rank*my_job;
  my_zone_end   = my_zone_start + my_job-1;
  // make sure last guy finishes it all
  if (my_rank == n_procs-1) my_zone_end = grid->z.size()-1;
  // get rid of unneeded processors
  if (my_zone_start >= grid->z.size()) my_zone_start = my_zone_end+1;

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
    if(verbose) cout << "ERROR: you must simulate at least one species of particle." << endl;
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

  // initialize all the zone eas variables
  #pragma omp parallel for collapse(2)
  for(int i=0; i<species_list.size(); i++) 
    for(int j=0; j<grid->z.size(); j++)
      species_list[i]->set_eas(j);

  // scatter initial particles in the simulation area
  // don't initialize them if iterative calculation. They all come from the core.
  int init_particles = lua->scalar<int>("init_particles");
  if(!iterate) initialize_particles(init_particles);


  //=================//
  // SET UP THE CORE //
  //=================//
  // the core temperature is used only in setting its emis vector
  // so it's looked at only in species::myInit()
  if(do_core){
    r_core      = lua->scalar<double>("r_core");
    L_core      = lua->scalar<double>("L_core");
    core_species_cdf.resize(species_list.size());
    for(int i=0; i<species_list.size(); i++){
      core_species_cdf.set_value(i, species_list[i]->int_core_emis());}
    core_species_cdf.normalize();
  }
  else{
    r_core = 0;
    L_core = 0;
  }

  // set the net luminosity
  L_net = L_core + therm_lum + decay_lum;
}




//------------------------------------------------------------
// take a transport time step 
//------------------------------------------------------------
void transport::step(double dt)
{
  // nominal time for iterative calc is 1
  if (iterate) dt = 1;
  
  // calculate the zone eas variables
  #pragma omp parallel for collapse(2)
  for(int i=0; i<species_list.size(); i++) 
    for(int j=0; j<grid->z.size(); j++)
      species_list[i]->set_eas(j);

  // emit new particles
  emit_particles(dt);

  // clear the tallies of the radiation quantities in each zone
  #pragma omp parallel for
  for (int i=0;i<grid->z.size();i++) 
  {
    grid->z[i].e_rad  = 0;
    grid->z[i].e_abs  = 0;
    grid->z[i].fx_rad = 0;
    grid->z[i].fy_rad = 0;
    grid->z[i].fz_rad = 0;
    grid->z[i].l_abs  = 0;
  }

  // Propagate the particles
  propagate_particles(dt);

  // properly normalize the radiative quantities
  #pragma omp parallel for
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

  // solve for T_gas and Ye structure if radiative eq. applied
  if (radiative_eq) solve_eq_zone_values();
   
  // advance time step
  if (!iterate) t_now += dt;
}




//----------------------------------------------------------------------------
// sum up the number of particles in all species
//----------------------------------------------------------------------------
int transport::total_particles(){
  return particles.size();
}


//----------------------------------------------------------------------------
// randomly sample the nu-integrated emissivities of all
// species to determine the species of a new particle
// emitted from the core
//----------------------------------------------------------------------------
int transport::sample_core_species()
{
  // randomly sample the species (precomputed)
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
int transport::sample_zone_species(int zone_index)
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
// get the doppler shift from lab to comoving
//------------------------------------------------------------
double transport::dshift_lab_to_comoving(particle* p)
{
  double v[3];
  grid->velocity_vector(p->ind,p->x,v);

  // outgoing velocity vector
  double beta = 0.0;
  for (int i=0;i<3;i++)  beta += v[i]*v[i];
  beta = sqrt(beta)/pc::c;
  double gamma  = 1.0/sqrt(1 - beta*beta);

  // doppler shifts outgoing
  double vdp    = (p->D[0]*v[0] + p->D[1]*v[1] + p->D[2]*v[2]);
  double dshift = gamma*(1 - vdp/pc::c);
  return dshift;
}


//------------------------------------------------------------
// get the doppler shift from comoving to lab
//------------------------------------------------------------
double transport::dshift_comoving_to_lab(particle* p)
{
  double v[3];
  grid->velocity_vector(p->ind,p->x,v);

  // outgoing velocity vector
  double beta = 0.0;
  for (int i=0;i<3;i++)
  {
    v[i] = -1*v[i];
    beta += v[i]*v[i];
  }
  beta = sqrt(beta)/pc::c;
  double gamma  = 1.0/sqrt(1 - beta*beta);

  // doppler shifts outgoing
  double vdp    = (p->D[0]*v[0] + p->D[1]*v[1] + p->D[2]*v[2]);
  double dshift = gamma*(1 - vdp/pc::c);
  return dshift;
}

//------------------------------------------------------------
// do a lorentz transformation; modifies the energy, frequency
// and direction vector of the particle
// sign =  1 for comoving to lab
// sign = -1 for lab to comoving
//------------------------------------------------------------

void transport::transform_comoving_to_lab(particle* p)
{
  lorentz_transform(p,1);
}
void transport::transform_lab_to_comoving(particle* p)
{
  lorentz_transform(p,-1);
}

void transport::lorentz_transform(particle* p, double sign)
{
  double v[3];
  grid->velocity_vector(p->ind,p->x,v);

  // outgoing velocity vector
  double beta = 0.0;
  for (int i=0;i<3;i++)
  {
    v[i] = sign*v[i];
    beta += v[i]*v[i];
  }
  beta = sqrt(beta)/pc::c;
  double gamma  = 1.0/sqrt(1 - beta*beta);

  // doppler shifts outgoing
  double vdp    = (p->D[0]*v[0] + p->D[1]*v[1] + p->D[2]*v[2]);
  double dshift = gamma*(1 - vdp/pc::c);

  // doppler shift the energy and frequency
  p->e   *= dshift;
  p->nu  *= dshift;

  // transform direction
  double D_old[3];
  D_old[0] = p->D[0];
  D_old[1] = p->D[1];
  D_old[2] = p->D[2];

  // See Mihalas & Mihalas eq 89.8
  p->D[0] = 1.0/dshift*(D_old[0]-gamma*v[0]/pc::c*(1-gamma*vdp/pc::c/(gamma+1)));
  p->D[1] = 1.0/dshift*(D_old[1]-gamma*v[1]/pc::c*(1-gamma*vdp/pc::c/(gamma+1)));
  p->D[2] = 1.0/dshift*(D_old[2]-gamma*v[2]/pc::c*(1-gamma*vdp/pc::c/(gamma+1)));

  // for security, make sure it is properly normalized
  double norm = p->D[0]*p->D[0] + p->D[1]*p->D[1] + p->D[2]*p->D[2];
  p->D[0] = p->D[0]/norm;
  p->D[1] = p->D[1]/norm;
  p->D[2] = p->D[2]/norm;
}
