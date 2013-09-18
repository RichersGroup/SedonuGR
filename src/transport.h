#ifndef _TRANSPORT_H
#define _TRANSPORT_H
#include "particle.h"
#include "spectrum_array.h"
#include "Lua.h"
#include "grid_general.h"
#include "locate_array.h"
#include "cdf_array.h"
#include "thread_RNG.h"
#include <list>


#define MAX_PARTICLES 1000000
class species_general;

class transport
{

private:

  // this species' list of particles
  std::vector<particle> particles;

  // the porition of zones this process is responsible for
  int my_zone_start, my_zone_end;

  // creation of particles functions
  void emit_particles(double dt);
  void emit_inner_source(double dt);
  void initialize_particles(int n_parts);
  void create_isotropic_particle(int zone_index, double Ep);
  int sample_core_species();
  int sample_zone_species(int zone_index);
  cdf_array core_species_cdf;

  // transformation functions
  void   lorentz_transform(particle* p, double);
  double dshift_comoving_to_lab(particle* p);
  double dshift_lab_to_comoving(particle* p);

  // propagate the particles
  void   transform_comoving_to_lab(particle* p);
  void   transform_lab_to_comoving(particle* p);

  // propagate the particles
  void propagate_particles(double dt);

  // scattering functions
  void propagate(particle* p, double dt);
  void isotropic_scatter(particle* p, int redistribute);

  // solve for temperature
  void   solve_eq_zone_values();
  double brent_method(int zone_index, double (*eq_function)(int,double,transport*), double min, double max);

public:

  // arrays of species
  vector<species_general*> species_list;

  // random number generator
  thread_RNG rangen;

  // current time in simulation
  double t_now;

  // remember what we're simulating
  int do_photons;
  int do_neutrinos;

  // remember what variables we're solving for
  int solve_T;
  int solve_Ye;

  // pointer to grid
  grid_general *grid;

  // items for core emission
  int n_inject;
  double r_core;
  double L_core;

  // simulation parameters
  double step_size;
  int    radiative_eq;
  int    iterate;
  int    verbose;
  double damping;

  // set things up
  void   init(Lua* lua);

  // in-simulation functions accessible to main
  void   step(double dt);
  int    total_particles();

  // stored minimum and maximum values for use by the Brent solver
  double T_min,  T_max;
  double Ye_min, Ye_max;
};

#endif

