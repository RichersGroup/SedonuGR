#ifndef _TRANSPORT_H
#define _TRANSPORT_H
#include "particle.h"
#include "spectrum_array.h"
#include "Lua.h"
#include "grid_general.h"
#include "locate_array.h"
#include "cdf_array.h"
#include <gsl/gsl_rng.h>
//#include "species_general.h"

#define MAX_PARTICLES 1000000
class species_general;

class transport
{

private:

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

  // solve for temperature
  void   solve_eq_zone_values();
  double brent_method(int zone_index, double (*eq_function)(int,double,transport*), double min, double max);
  //double temp_eq_function(int zone_index, double T);
  //double Ye_eq_function(int zone_index, double T);

public:

  // arrays of species
  vector<species_general*> species_list;

  // random number generator
  gsl_rng *rangen;

  // current time in simulation
  double t_now;

  // remember what we're simulating
  int do_photons;
  int do_neutrinos;

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

  // set things up
  void   init(Lua* lua);

  // in-simulation functions accessible to main
  void   step(double dt);
  int    total_particles();
  void   update_composition();

  // stored minimum and maximum values for use by the Brent solver
  double T_min,  T_max;
  double Ye_min, Ye_max;
};

#endif

