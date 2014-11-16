#ifndef _TRANSPORT_H
#define _TRANSPORT_H
#include <mpi.h>
#include <vector>
#include "particle.h"
#include "Lua.h"
#include "locate_array.h"
#include "cdf_array.h"
#include "thread_RNG.h"
#include "global_options.h"

class species_general;
class grid_general;
class spectrum_array;
enum ParticleEvent {interact, zoneEdge, timeStep, boundary};

class transport
{

private:

	// this species' list of particles
	std::vector<particle> particles;

	// MPI stuff
	int MPI_nprocs;
	int MPI_myID;
	void reduce_radiation();
	void synchronize_gas();
	vector<unsigned> my_zone_end;

	// subroutine for calculating timescales
	void calculate_timescales() const;
	void calculate_annihilation() const;

	// main function to emit particles
	void emit_particles(double dt);

	// emit from where?
	//void initialize_particles(int init_particles);
	void initialize_blackbody(const double T, const double munue);
	void emit_inner_source(const int n_emit, const double dt, double t=-1);
	void emit_zones(const int n_emit, const double dt, double t=-1);

	// what kind of particle to create?
	void create_surface_particle(const double Ep, const double t);
	void create_thermal_particle(const int zone_index, const double Ep, const double t);

	// per-zone luminosity functions
	double zone_comoving_visc_heat_rate(const int zone_index) const;
	double zone_comoving_therm_emit_energy (const int zone_index, const double lab_dt) const;
	double zone_comoving_therm_emit_leptons(const int zone_index, const double lab_dt) const;

	// species sampling functions
	int sample_core_species() const;
	int sample_zone_species(int zone_index) const;

	// transformation functions
	double comoving_dt(const double lab_dt, const int z_ind) const;
	double dshift_comoving_to_lab   (const particle* p, const int z_ind=-1) const;
	double dshift_lab_to_comoving   (const particle* p, const int z_ind=-1) const;
	void   transform_comoving_to_lab(particle* p, const int z_ind=-1) const;
	void   transform_lab_to_comoving(particle* p, const int z_ind=-1) const;

	// propagate the particles
	void propagate_particles(const double dt);
	void propagate(particle* p, const double dt);
	void move(particle* p, const double lab_d);
	void tally_radiation(const particle* p, const int z_ind, const double dshift_l2c, const double lab_d, const double lab_opac, const double abs_frac) const;
	void reset_radiation();
	void which_event(const particle* p,const int z_ind, const double dt, const double lab_opac, double* d_smallest, ParticleEvent *event) const;
	void event_boundary(particle* p, const int z_ind) const;
	void event_interact(particle* p, const int z_ind, const double abs_frac);
	void isotropic_scatter(particle* p) const;
	void re_emit(particle* p, const int z_ind) const;
	void remove_dead_particles();

	// solve for temperature and Ye (if steady_state)
	int    solve_T;
	int    solve_Ye;
	double damping;
	int    brent_itmax;
	double brent_solve_tolerance;
	void   solve_eq_zone_values();
	void   normalize_radiative_quantities(const double dt);
	double brent_method(const int zone_index, double (*eq_function)(int,double,transport*), const double min, const double max);

	// update temperature and Ye (if !steady_state)
	void update_zone_quantities();

	// stored minimum and maximum values to assure safety
	double T_min,  T_max;
	double Ye_min, Ye_max;
	double rho_min, rho_max;
	int max_particles;

	// simulation parameters
	double step_size;
	int    do_photons;
	int    do_neutrinos;
	int    do_distribution;
	int    iterative;
	int    radiative_eq;
	int    rank0;
	double ratio_emit_by_zone;

	// output parameters
	int write_zones_every;
	int write_rays_every;
	int write_spectra_every;

	// current time in simulation
	double t_now;

	// check parameters
	void check_parameters() const;

	// check quality
	double Q_zones() const;
	double Q_core() const;

public:

	transport();

	int verbose;
	double current_time();

	// arrays of species
	vector<species_general*> species_list;

	// pointer to grid
	grid_general *grid;

	// items for core emission
	double r_core;
	int n_emit_core;
	double core_lum_multiplier;
	cdf_array core_species_luminosity;
	void init_core(const double r_core, const double T_core, const double munue_core);

	// items for zone emission
	int do_visc;
	int n_emit_zones;
	double visc_specific_heat_rate;

	// initial particle creation
	int n_initial;
	double initial_BB_T;
	double initial_BB_munue;

	// how many times do we emit+propagate each timestep?
	int emissions_per_timestep;

	// global radiation quantities
	double L_net_lab;
	double L_esc_lab;
	vector<long> n_active;
	vector<long> n_escape;


	// reflect off the outer boundary?
	int reflect_outer;

	// random number generator
	mutable thread_RNG rangen;

	// set things up
	void   init(Lua* lua);

	// in-simulation functions to be used by main
	void step(const double dt);
	void write(const int it) const;
	int  total_particles() const;
	void write_rays(const int it);
	static void open_file(const char* filebase, const int iw, ofstream& outf);
	static double lorentz_factor(const vector<double>& v);
	static double dot(const vector<double>& a, const vector<double>& b);
	static void normalize(vector<double>& a);
	static double mean_mass(const double Ye);
};

#endif

