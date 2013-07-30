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
#include "transport.h"
#include "Lua.h"
//#include "radioactive.hh"
//#include <omp.h>

namespace pc = physical_constants;



//--------------------------------------------------------
// constructor
//--------------------------------------------------------
transport::transport()
{
  // defaults
  step_size = 0.1;
  grey_opac = 0.0;
  r_core = 0;
}


//--------------------------------------------------------
// Initialize and allocate
//--------------------------------------------------------
void transport::init(string infile,grid_general *g)
{
  // save param file name
  param_file = infile;
  

  // start at time 0
  t_now = 0;

  // point to the grid
  this->grid = g;
  
  // get mpi rank, verbosity
  int my_rank, n_procs;
  MPI_Comm_size( MPI_COMM_WORLD, &n_procs);
  MPI_Comm_rank( MPI_COMM_WORLD, &my_rank );
  if (my_rank == 0) verbose = 1; else verbose = 0;

  // figure out which zones are in this processors work load
  int my_job = (int)(grid->n_zones/(1.0*n_procs));
  if (my_job < 1) my_job = 1;
  my_zone_start = my_rank*my_job;
  my_zone_end   = my_zone_start + my_job;
  // make sure last guy finishes it all
  if (my_rank == n_procs-1) my_zone_end = grid->n_zones;
  // get rid of unneeded processers
  if (my_zone_start >= grid->n_zones) my_zone_start = my_zone_end;


  // setup and seed random number generator
  const gsl_rng_type * TypeR;
  gsl_rng_env_setup();
  gsl_rng_default_seed = (unsigned int)time(NULL) + my_rank;
  TypeR = gsl_rng_default;
  rangen = gsl_rng_alloc (TypeR);

  // open up the lua parameter file
  Lua lua;
  lua.init(infile);

  // read relevant parameters
  this->radiative_eq = lua.scalar<int>("radiative_eq");
  this->iterate      = lua.scalar<int>("iterate");
  this->step_size = lua.scalar<double>("step_size");
  this->epsilon   = lua.scalar<double>("epsilon");

  
  // intialize output spectrum
  std::vector<double>stg = lua.vector<double>("spec_time_grid");
  std::vector<double>sng = lua.vector<double>("spec_nu_grid");
  int nmu  = lua.scalar<int>("n_mu");
  int nphi = lua.scalar<int>("n_phi");
  spectrum.init(stg,sng,nmu,nphi);
  spectrum.set_name("optical_spectrum.dat");

  // initalize opacities
  initialize_opacity(&lua);
  set_opacity();
  
  // initialize particles
  int n_parts = lua.scalar<int>("init_particles");
  initialize_particles(n_parts);

  // close lua file
  lua.close();
}




//------------------------------------------------------------
// take a transport time step 
//------------------------------------------------------------
void transport::step(double dt)
{
  int i;

  // nominal time for iterative calc is 1
  if (iterate) dt = 1;
  
  // calculate opacities
  set_opacity();

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
  vector<particle>::iterator pIter = particles.begin();
  int n_active = particles.size();
  int n_escape = 0;
  while (pIter != particles.end())
  {
    ParticleFate fate = propagate(*pIter,dt);
    if (fate == escaped) n_escape++;
    if ((fate == escaped)||(fate == absorbed)) particles.erase(pIter);
    else pIter++;
  }
  double per_esc = (100.0*n_escape)/n_active;
  if ((verbose)&&(iterate)) cout << "# Percent escaped = " << per_esc << endl;
  spectrum.rescale(1.0/per_esc);
  
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
    double d_bn = step_size*grid->zone_min_length(p.ind);

    // doppler shift from comoving to lab
    double dshift = dshift_comoving_to_lab(p);

    // get local opacity and absorption fraction (epsilon)
    double opac, eps;
    this->get_opacity(p,dshift,opac,eps);
    
    // convert opacity from comoving to lab frame for the purposes of 
    // determining the interaction distance in the lab frame
    // This corresponds to equation 90.8 in Mihalas&Mihalas. You multiply 
    // the comoving opacity by nu_0 over nu, which is why you
    // multiply by dshift instead of dividing by dshift here
    opac = opac*dshift;

    // random optical depth to next interaction
    tau_r = -1.0*log(1 - gsl_rng_uniform(rangen));
    
    // step size to next interaction event
    d_sc  = tau_r/opac;
    if (opac == 0) d_sc = INFINITY;
    if (d_sc < 0) cout << "ERROR: negative interaction distance!\n";
  
    // find distance to end of time step
    d_tm = (tstop - p.t)*pc::c;
    // if iterative calculation, let all particles escape
    if (iterate) d_tm = INFINITY;

    // find out what event happens (shortest distance)
    if ((d_sc < d_bn)&&(d_sc < d_tm))
      {event = scatter;    this_d = d_sc;}
    else if (d_bn < d_tm)
      {event = boundary;   this_d = d_bn;}
    else 
      {event = tstep;      this_d = d_tm; }

    // tally in contribution to zone's radiation energy (both *lab* frame)
    double this_E = p.e*this_d; 
    zone->e_rad += this_E; 

    // shift opacity back to comoving frame for energy and momentum exchange. 
    // Radiation energy is still lab frame
    opac = opac / dshift;

    // store absorbed energy in *comoving* frame 
    // (will turn into rate by dividing by dt later)
    // Extra dshift definitely needed here (two total)
    zone->e_abs  += this_E*dshift*(opac)*eps*zone->eps_imc * dshift; 

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
      double z = gsl_rng_uniform(rangen);
      
      // do photon interaction physics
      if (p.type == photon)
      {
	// see if scattered 
	if (z > eps) isotropic_scatter(p,0);
	else
	{
	  // check for effective scattering
	  double z2;
	  if (radiative_eq) z2 = 2;
	  else z2 = gsl_rng_uniform(rangen);
	  // do an effective scatter
	  if (z2 > zone->eps_imc) isotropic_scatter(p,1);
	  // otherwise really absorb (kill) it
	  else fate = absorbed; 
	}
      }
    }
    
    // ---------------------------------
    // do if time step end
    // ---------------------------------
    else if (event == tstep) { fate = stopped;}
    
    // Find position of the particle now
    int ix[3];
    p.ind = grid->get_zone(p.x);
    if (p.ind == -1) fate = absorbed;
    if (p.ind == -2) fate = escaped;

    // check for inner boundary absorption
    if (p.r() < r_core)  {fate = absorbed;}
  }

  // Add escaped photons to output spectrum
  if (fate == escaped) 
    if (p.type == photon)
    {
      // account for light crossing time, relative to grid center
      //double X0 = p.x[0] - grid->x_cen;
      //double X1 = p.x[1] - grid->x_cen;
      //double X2 = p.x[2] - grid->x_cen;
      
      //double xdot = X0*p.D[0] + X1*p.D[1] + X2*p.D[2];
      double t_obs = p.t; // - xdot/pc::c;
      spectrum.count(t_obs,p.nu,p.e,p.D);
    }

  return fate;
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


//***************************************************************/
// This is the function that expresses radiative equillibrium
// in a cell (i.e. E_absorbed = E_emitted).  It is used in
// The Brent solver below to determine the temperature such
// that RadEq holds
//**************************************************************/
double transport::rad_eq_function(int c,double T)
{
  // total energy absorbed in zone
  double E_absorbed = grid->z[c].e_abs;
  // total energy emitted (to be calculated)
  double E_emitted = 0.;

  // integrate emisison over frequency (angle
  // integration gives the 4*PI) to get total
  // radiation energy emitted. Opacities are
  // held constant for this (assumed not to change
  // much from the last time step).
  for (int i=0;i<nu_grid.size();i++)
  {
    double dnu  = nu_grid.delta(i);
    double nu   = nu_grid.center(i);
    double B_nu = blackbody_nu(T,nu);
    double kappa_abs  = epsilon*grid->z[c].opac[i];
    E_emitted += 4.0*pc::pi*kappa_abs*B_nu*dnu;
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


