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

class species_general;

class transport
{

private:

  // this species' list of particles
  std::vector<particle> particles;

  // the porition of zones this MPI process is responsible for
  int my_zone_start, my_zone_end;

  // creation of particles functions
  void emit_particles(double dt);
  void emit_inner_source(double dt);
  void emit_zones(double dt);
  void initialize_particles(int init_particles);
  void create_surface_particle(double Ep, double t);
  void create_thermal_particle(int zone_index, double Ep, double t);
  void create_decay_particle(int zone_index, double Ep, double t);
  int sample_core_species();
  int sample_zone_species(int zone_index);
  double zone_visc_heat_rate(int zone_index);
  double zone_heat_lum(int zone_index);
  double zone_decay_lum(int zone_index);

  // transformation functions
  void   lorentz_transform        (particle* p, double);
  double dshift_comoving_to_lab   (particle* p);
  double dshift_lab_to_comoving   (particle* p);
  void   transform_comoving_to_lab(particle* p);
  void   transform_lab_to_comoving(particle* p);

  // propagate the particles
  void propagate_particles(double dt);
  void propagate(particle* p, double dt);
  void isotropic_scatter(particle* p, int redistribute);

  // solve for temperature and Ye
  void   solve_eq_zone_values();
  double brent_method(int zone_index, double (*eq_function)(int,double,transport*), double min, double max);

  // stored minimum and maximum values to assure safety
  double T_min,  T_max;
  double Ye_min, Ye_max;
  int max_particles;

  // simulation parameters
  double step_size;
  double damping;
  int    do_photons;
  int    do_neutrinos;
  int    solve_T;
  int    solve_Ye;
  int    radiative_eq;
  int    iterate;
  int    verbose;

  // current time in simulation
  double t_now;

public:

  // arrays of species
  vector<species_general*> species_list;

  // pointer to grid
  grid_general *grid;

  // items for core emission
  bool do_core;
  double r_core;
  double L_core;
  int n_emit_core;
  cdf_array core_species_cdf;

  // items for zone emission
  bool do_therm, do_decay;
  int n_emit_therm, n_emit_decay;
  double visc_specific_heat_rate;

  // net luminosity
  double L_net;

  // random number generator
  thread_RNG rangen;

  // set things up
  void   init(Lua* lua);

  // in-simulation functions to be used by main
  void   step(double dt);
  int    total_particles();
};

#endif

