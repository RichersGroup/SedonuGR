#ifndef _TRANSPORT_H
#define _TRANSPORT_H
#include "particle.h"
#include "spectrum_array.h"
#include "Lua.h"
#include "grid_general.h"
#include "locate_array.h"
#include "cdf_array.h"

#define MAX_PARTICLES 1000000
#define TEMP_RANGE_MAX 1.e12
#define TEMP_RANGE_MIN 1.

enum ParticleFate  {moving, stopped, escaped, absorbed};

class transport
{

private:

  // arrays of particles
  vector<particle> particles;

   // random number generator
  gsl_rng *rangen;

  // name of the input parameter file
  string param_file;

  // the porition of zones this process is responsible for
  int my_zone_start, my_zone_end;

  // opacity functions
  void   initialize_opacity(Lua*);
  void   get_opacity(particle &p, double dshift, double &opac, double &eps);
  void   set_opacity();
  double klein_nishina(double);
  double blackbody_nu(double T, double nu);

  // propagate the particles functions
  ParticleFate propagate(particle &p, double tstop);
  void   lorentz_transform(particle &p, double);
  void   transform_comoving_to_lab(particle &p);
  void   transform_lab_to_comoving(particle &p);
  double dshift_comoving_to_lab(particle p);
  double dshift_lab_to_comoving(particle p);

  // creation of particles functions
  void   emit_particles(double dt);
  void   emit_inner_source(double dt);
  void   create_isotropic_particle(int,double);

  // scattering functions
  void   compton_scatter(particle*);
  void   isotropic_scatter(particle &p, int);

  // solve for temperature
  void   solve_eq_temperature();
  double temp_brent_method(int zone);
  double rad_eq_function(int c, double T);

  // frequency grid, describing the frequency spacing
  locate_array nu_grid;
  
  // items for core emission
  int n_inject;
  double r_core;
  cdf_array core_emis;


public:

  // current time in simulation
  double t_now;
  
  // pointer to grid
  grid_general *grid;

  // emmergent spectrum,
  spectrum_array spectrum;

  int    verbose;
  double n_photons_per;         // number of photons to emit per day
  double grey_opac;
  double emit_min, emit_max;
  double step_size;
  double epsilon;
  int radiative_eq;
  int iterate;
 
  transport();
  void   init(string infile, grid_general*);
  void   initialize_particles(int);

  void   step(double dt);
  int    num_particles()        {return particles.size();}
 
};

#endif

