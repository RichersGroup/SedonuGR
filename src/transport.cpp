#include <omp.h>
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
    if(verbose) std::cout << "ERROR: the requested grid type is not implemented." << std::endl;
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
  my_zone_end   = my_zone_start + my_job-1;
  // make sure last guy finishes it all
  if (my_rank == n_procs-1) my_zone_end = grid->z.size()-1;
  // get rid of unneeded processors
  if (my_zone_start >= grid->z.size()) my_zone_start = my_zone_end+1;

  // setup and seed random number generator
  rangen.init();

  //==================//
  // SET UP TRANSPORT //
  //==================//
  // start at time 0
  t_now = 0;

  // read relevant parameters
  radiative_eq  = lua->scalar<int>("radiative_eq");
  iterate       = lua->scalar<int>("iterate");
  step_size     = lua->scalar<double>("step_size");
  damping       = lua->scalar<double>("damping");
  solve_T       = lua->scalar<int>("solve_T");
  solve_Ye      = lua->scalar<int>("solve_Ye");


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
  // solve radiative equilibrium temperature and Ye (but only for the zones I'm responsible for)
  #pragma omp parallel for schedule(guided)
  for (int i=0; i<grid->z.size(); i++)
  {
    // set Ye and temp to zero so it's easy to reduce later
    if(i<my_zone_start || i>my_zone_end){
      grid->z[i].T_gas = 0;
      grid->z[i].Ye    = 0;
      continue;  // don't actually solve for these zones
    }

    double T_last_iter, Ye_last_iter;
    double T_last_step, Ye_last_step;
    double T_error,Ye_error;
    double dT_step, dYe_step;
    int iter=0;

    iter = 0;

    // set up the solver
    if(solve_T)
    {
      T_error  = 10*BRENT_SOLVE_TOLERANCE;
      T_last_step  = grid->z[i].T_gas;
    }
    if(solve_Ye)
    {
      Ye_error = 10*BRENT_SOLVE_TOLERANCE;
      Ye_last_step = grid->z[i].Ye;
    }

    // loop through solving the temperature and Ye until both are within error.
    while(iter<=BRENT_ITMAX && (T_error>BRENT_SOLVE_TOLERANCE || Ye_error>BRENT_SOLVE_TOLERANCE))
    {
      if(solve_T)
      {
	T_last_iter  = grid->z[i].T_gas;
	grid->z[i].T_gas = brent_method(i, temp_eq_function, T_min,  T_max);
	T_error  = fabs( (grid->z[i].T_gas - T_last_iter ) / (T_last_iter ) );
      }
      if(solve_Ye)
      {
	Ye_last_iter = grid->z[i].Ye;
	grid->z[i].Ye    = brent_method(i, Ye_eq_function,   Ye_min, Ye_max);
	Ye_error = fabs( (grid->z[i].Ye    - Ye_last_iter) / (Ye_last_iter) );
      }
      iter++;
    }

    // warn if it didn't converge
    if(iter == BRENT_ITMAX) cout << "# WARNING: outer Brent solver hit maximum iterations." << endl;

    // damp the oscillations between steps, ensure that it's within the allowed boundaries
    if(solve_T)
    {
      dT_step  = grid->z[i].T_gas - T_last_step;
      grid->z[i].T_gas =  T_last_step + (1.0 - damping)*dT_step;
      if(grid->z[i].T_gas > T_max){
	cout << "# WARNING: Changing T_gas in zone " << i << " from " << grid->z[i].T_gas << " to T_max=" << T_max << endl;
	grid->z[i].T_gas = T_max;}
      if(grid->z[i].T_gas < T_min){
	cout << "# WARNING: Changing T_gas in zone " << i << " from " << grid->z[i].T_gas << " to T_min=" << T_min << endl;
	grid->z[i].T_gas = T_min;}
    }
    if(solve_Ye)
    {
      dYe_step = grid->z[i].Ye - Ye_last_step;
      grid->z[i].Ye = Ye_last_step + (1.0 - damping)*dYe_step;
      if(grid->z[i].Ye > Ye_max){
	cout << " WARNING: Changing Ye in zone " << i << " from " << grid->z[i].Ye << " to Ye_max=" << Ye_max << endl;
	grid->z[i].Ye = Ye_max;}
      if(grid->z[i].Ye < Ye_min){
	cout << " WARNING: Changing Ye in zone " << i << " from " << grid->z[i].Ye << " to Ye_min=" << Ye_min << endl;
	grid->z[i].Ye = Ye_min;}
    }
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
double transport::dshift_lab_to_comoving(particle p)
{
  double v[3];
  grid->velocity_vector(p.ind,p.x,v);

  // outgoing velocity vector
  double beta = 0.0;
  for (int i=0;i<3;i++)  beta += v[i]*v[i];
  beta = sqrt(beta)/pc::c;
  double gamma  = 1.0/sqrt(1 - beta*beta);

  // doppler shifts outgoing
  double vdp    = (p.D[0]*v[0] + p.D[1]*v[1] + p.D[2]*v[2]);
  double dshift = gamma*(1 - vdp/pc::c);
  return dshift;
}


//------------------------------------------------------------
// get the doppler shift from comoving to lab
//------------------------------------------------------------
double transport::dshift_comoving_to_lab(particle p)
{
  double v[3];
  grid->velocity_vector(p.ind,p.x,v);

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
  double vdp    = (p.D[0]*v[0] + p.D[1]*v[1] + p.D[2]*v[2]);
  double dshift = gamma*(1 - vdp/pc::c);
  return dshift;
}

//------------------------------------------------------------
// do a lorentz transformation; modifies the energy, frequency
// and direction vector of the particle
// sign =  1 for comoving to lab
// sign = -1 for lab to comoving
//------------------------------------------------------------

void transport::transform_comoving_to_lab(particle &p)
{
  lorentz_transform(p,1);
}
void transport::transform_lab_to_comoving(particle &p)
{
  lorentz_transform(p,-1);
}

void transport::lorentz_transform(particle &p, double sign)
{
  double v[3];
  grid->velocity_vector(p.ind,p.x,v);

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
  double vdp    = (p.D[0]*v[0] + p.D[1]*v[1] + p.D[2]*v[2]);
  double dshift = gamma*(1 - vdp/pc::c);

  // doppler shift the energy and frequency
  p.e   *= dshift;
  p.nu  *= dshift;

  // transform direction
  double D_old[3];
  D_old[0] = p.D[0];
  D_old[1] = p.D[1];
  D_old[2] = p.D[2];

  // See Mihalas & Mihalas eq 89.8
  p.D[0] = 1.0/dshift*(D_old[0]-gamma*v[0]/pc::c*(1-gamma*vdp/pc::c/(gamma+1)));
  p.D[1] = 1.0/dshift*(D_old[1]-gamma*v[1]/pc::c*(1-gamma*vdp/pc::c/(gamma+1)));
  p.D[2] = 1.0/dshift*(D_old[2]-gamma*v[2]/pc::c*(1-gamma*vdp/pc::c/(gamma+1)));

  // for security, make sure it is properly normalized
  double norm = p.D[0]*p.D[0] + p.D[1]*p.D[1] + p.D[2]*p.D[2];
  p.D[0] = p.D[0]/norm;
  p.D[1] = p.D[1]/norm;
  p.D[2] = p.D[2]/norm;
}


void transport::propagate_particles(double dt)
{
  // OPTIMIZE - implement forward list. stores 1 pointer instead of 2 in each element
  // TODO - the parallel algorithm causes lots of contention over the RNG. Need a separate RNG for each thread
  vector<long> n_active(species_list.size(),0);
  vector<long> n_escape(species_list.size(),0);
  double e_esc = 0;

  #pragma omp parallel default(none) shared(n_active,n_escape,dt,e_esc)
  #pragma omp single
  {
    list<particle>::iterator tmpIter;
    list<particle>::iterator pIter = particles.begin();

    // create a task to move each particle
    while(pIter != particles.end()){
      tmpIter = pIter;  // done like this to prevent race condition when incrementing pIter
      pIter++;
      n_active[tmpIter->s]++;

      //================================================
      #pragma omp task default(none) firstprivate(tmpIter) shared(dt,n_escape,e_esc)
      {
	ParticleFate fate = propagate(*tmpIter,dt);
	if (fate == escaped){ 
          #pragma omp atomic
	  n_escape[tmpIter->s]++;
          #pragma omp atomic
	  e_esc += tmpIter->e;
	  species_list[tmpIter->s]->spectrum.count(tmpIter->t, tmpIter->nu, tmpIter->e, tmpIter->D);
	}
	if ((fate == escaped)||(fate == absorbed)){
          #pragma omp critical
	  particles.erase(tmpIter);
	}
      }
      //================================================

    } //while
  } //#pragma omp parallel
  
  for(int i=0; i<species_list.size(); i++){
    double per_esc = (100.0*n_escape[i])/n_active[i];
    if (verbose && iterate){
      if(n_active[i]>0) cout << "# " << per_esc << "% of the " << species_list[i]->name << " escaped (" 
			     << n_escape[i] << "/" << n_active[i] << ")." << endl;
      else cout << "# " << "No active " << species_list[i]->name << endl;
    }
    if(per_esc>0) species_list[i]->spectrum.rescale(100.0/per_esc);
  }
  if(verbose && iterate) cout << "# Energy escaped: " << e_esc << " ergs." << endl;
}

//--------------------------------------------------------
// Propagate a single monte carlo particle until
// it  escapes, is absorbed, or the time step ends
//--------------------------------------------------------
ParticleFate transport::propagate(particle &p, double dt)
{
  enum ParticleEvent {scatter, boundary, tstep};
  ParticleEvent event;

  ParticleFate  fate = moving;

  // time of end of timestep
  double tstop = t_now + dt;

  // local variables
  double tau_r,d_sc,d_tm,this_d;

  // pointer to current zone
  zone *zone = &(grid->z[p.ind]);

  // propagate until this flag is set
  while (fate == moving)
  {
    // set pointer to current zone
    zone = &(grid->z[p.ind]);

    // maximum step size inside zone
    double d_bn = step_size * grid->zone_min_length(p.ind);

    // doppler shift from comoving to lab
    double dshift = dshift_comoving_to_lab(p);

    // get local opacity and absorption fraction
    double opac, abs_frac;
    species_list[p.s]->get_opacity(p,dshift,&opac,&abs_frac);

    // convert opacity from comoving to lab frame for the purposes of
    // determining the interaction distance in the lab frame
    // This corresponds to equation 90.8 in Mihalas&Mihalas. You multiply
    // the comoving opacity by nu_0 over nu, which is why you
    // multiply by dshift instead of dividing by dshift here
    double opac_lab = opac*dshift;

    // random optical depth to next interaction
    tau_r = -1.0*log(1 - rangen.uniform());

    // step size to next interaction event
    d_sc  = tau_r/opac_lab;
    if (opac_lab == 0) d_sc = INFINITY;
    if (d_sc < 0){
      cout << "ERROR: negative interaction distance!\n" << endl;
      cout << __FILE__ << ":" << __LINE__ << endl;
      exit(15);
    }

    // find distance to end of time step
    d_tm = (tstop - p.t)*pc::c;
    if (iterate) d_tm = INFINITY; // i.e. let all particles escape

    // find out what event happens (shortest distance)
    if ( (d_sc < d_bn) && (d_sc < d_tm) ){
      event  = scatter;
      this_d = d_sc;
    }
    else if (d_bn < d_tm){
      event  = boundary;
      this_d = d_bn;
    }
    else{
      event  = tstep;
      this_d = d_tm;
    }

    // tally in contribution to zone's radiation energy (both *lab* frame)
    double this_E = p.e*this_d;
    #pragma omp atomic
    zone->e_rad += this_E;

    // store absorbed energy in *comoving* frame
    // (will turn into rate by dividing by dt later)
    // Extra dshift definitely needed here (two total)
    // to convert both p.e and this_d to the comoving frame
    double this_E_comoving = this_E * dshift * dshift;
    #pragma omp atomic
    zone->e_abs += this_E_comoving * (opac*abs_frac*zone->eps_imc);

    // store absorbed lepton number (same in both frames, except for the
    // factor of this_d which is divided out later
    if(species_list[p.s]->lepton_number != 0){
      double this_l_comoving = species_list[p.s]->lepton_number * p.e/(p.nu*pc::h) * this_d*dshift;
      #pragma omp atomic
      zone->l_abs += this_l_comoving * (opac*abs_frac*zone->eps_imc);
    }

    // put back in radiation force tally here
    // fx_rad =

    // move particle the distance
    p.x[0] += this_d*p.D[0];
    p.x[1] += this_d*p.D[1];
    p.x[2] += this_d*p.D[2];
    // advance the time
    p.t = p.t + this_d/pc::c;

    // ---------------------------------
    // Do if scatter
    // ---------------------------------
    if (event == scatter)
    {
      // random number to check for scattering or absorption
      double z = rangen.uniform();

      // decide whether to scatter
      if (z > abs_frac) isotropic_scatter(p,0);
      // or absorb
      else
      {
	// check for effective scattering
	double z2;
	if (radiative_eq) z2 = 2;
	else z2 = rangen.uniform();

	// do an effective scatter (i.e. particle is absorbed
	// but, since we require energy in = energy out it is re-emitted)
	if (z2 > zone->eps_imc) isotropic_scatter(p,1);
	// otherwise really absorb (kill) it
	else fate = absorbed;
      }
    }

    // ---------------------------------
    // do if time step end
    // ---------------------------------
    else if (event == tstep) fate = stopped;

    // Find position of the particle now
    p.ind = grid->get_zone(p.x);
    if (p.ind == -1) fate = absorbed;
    if (p.ind == -2) fate = escaped;

    // check for inner boundary absorption
    if (p.r() < r_core) fate = absorbed;
  }

  return fate;
}


