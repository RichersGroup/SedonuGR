#ifndef _TRANSPORT_H
#define _TRANSPORT_H
#include <mpi.h>
#include <vector>
#include "particle.h"
#include "spectrum_array.h"
#include "Lua.h"
#include "grid_general.h"
#include "locate_array.h"
#include "cdf_array.h"
#include "thread_RNG.h"

class species_general;

class transport
{

private:

  // this species' list of particles
  std::vector<particle> particles;

  // MPI stuff
  int MPI_nprocs;
  int MPI_myID;
  MPI_Datatype MPI_real;
  void reduce_radiation();
  void synchronize_gas();
  vector<int> my_zone_end;

  // main function to emit particles
  void emit_particles(double dt);

  // emit from where?
  //void initialize_particles(int init_particles);
  void emit_inner_source(double dt);
  void emit_zones(const double dt, 
		  const int n_emit, 
		  double (transport::*zone_lum)(const int) const, 
		  void (transport::*emit_particle)(const int,const double,const double));

  // what kind of particle to create?
  void create_surface_particle(const double Ep, const double t);
  void create_thermal_particle(const int zone_index, const double Ep, const double t);
  void create_decay_particle(const int zone_index, const double Ep, const double t);

  // per-zone luminosity functions
  double zone_visc_heat_rate(const int zone_index) const;
  double zone_therm_lum(const int zone_index) const;
  double zone_decay_lum(const int zone_index) const;
  
  // species sampling functions
  int sample_core_species() const;
  int sample_zone_species(int zone_index) const;

  // transformation functions
  double dshift_comoving_to_lab   (particle* p) const;
  double dshift_lab_to_comoving   (particle* p) const;
  void   transform_comoving_to_lab(particle* p) const;
  void   transform_lab_to_comoving(particle* p) const;

  // propagate the particles
  void propagate_particles(const double dt);
  void propagate(particle* p, const double dt) const;
  void isotropic_scatter(particle* p, const int redistribute) const;

  // solve for temperature and Ye
  int    solve_T;
  int    solve_Ye;
  double damping;
  int    brent_itmax;
  double brent_solve_tolerance;
  void   solve_eq_zone_values();
  double brent_method(const int zone_index, double (*eq_function)(int,double,transport*), const double min, const double max);

  // stored minimum and maximum values to assure safety
  double T_min,  T_max;
  double Ye_min, Ye_max;
  double rho_min, rho_max;
  int max_particles;

  // simulation parameters
  double step_size;
  int    do_photons;
  int    do_neutrinos;
  int    steady_state;
  int    verbose;

  // current time in simulation
  double t_now;

public:

  // arrays of species
  vector<species_general*> species_list;

  // pointer to grid
  grid_general *grid;

  // items for core emission
  double r_core;
  double L_core;
  int n_emit_core;
  cdf_array core_species_cdf;

  // items for zone emission
  int do_visc;
  int n_emit_therm, n_emit_decay, n_emit_visc;
  double visc_specific_heat_rate;

  // net luminosity
  double L_net;

  // random number generator
  mutable thread_RNG rangen;

  // set things up
  void   init(Lua* lua);

  // in-simulation functions to be used by main
  void step(const double dt);
  int  total_particles() const;
  void write_zones(const int it);
  void write_rays(const int it);
  void write_spectra(const int it);
};

#endif

