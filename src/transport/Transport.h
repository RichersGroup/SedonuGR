/*
//  Copyright (c) 2015, California Institute of Technology and the Regents
//  of the University of California, based on research sponsored by the
//  United States Department of Energy. All rights reserved.
//
//  This file is part of Sedonu.
//
//  Sedonu is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  Neither the name of the California Institute of Technology (Caltech)
//  nor the University of California nor the names of its contributors 
//  may be used to endorse or promote products derived from this software
//  without specific prior written permission.
//
//  Sedonu is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with Sedonu.  If not, see <http://www.gnu.org/licenses/>.
//
*/

#ifndef _TRANSPORT_H
#define _TRANSPORT_H
#include <vector>
#include "Particle.h"
#include "Lua.h"
#include "CDFArray.h"
#include "ThreadRNG.h"
#include "LorentzHelper.h"

class Species;
class Grid;
class SpectrumArray;
enum ParticleEvent {interact, nothing};

class Transport
{

protected:

	// this species' list of particles
	std::vector<Particle> particles;

	// MPI stuff
	int MPI_nprocs;
	int MPI_myID;
	void reduce_radiation();
	void synchronize_gas();
	std::vector<unsigned> my_zone_end;

	// subroutine for calculating timescales
	void calculate_annihilation() const;

	// main function to emit particles
	void emit_particles();

	// emit from where?
	void emit_inner_source();
	void emit_zones();
	void emit_inner_source_by_bin();
	void emit_zones_by_bin();

	// what kind of particle to create?
	void create_surface_particle(const double Ep, const int s=-1, const int g=-1);
	void create_thermal_particle(const int zone_index, const double Ep, const int s=-1, const int g=-1);

	// species sampling functions
	int sample_core_species() const;
	int sample_zone_species(const int zone_index, double *Ep) const;

	// transformation functions
	void get_opacity(LorentzHelper *lh, const int z_ind) const;

	// propagate the particles
	void propagate_particles();
	void propagate(Particle* p);
	virtual void move(LorentzHelper *lh, int *z_ind);
	void tally_radiation(const LorentzHelper *lh, const int z_ind) const;
	void reset_radiation();
	void which_event(LorentzHelper* lh,const int z_ind, ParticleEvent *event) const;
	void boundary_conditions(LorentzHelper *lh, int *z_ind) const;
	void event_interact(LorentzHelper* lh, int *z_ind);
	void scatter(LorentzHelper *lh, int *z_ind) const;
	void random_walk(LorentzHelper *lh, const double Rcom, const double D, const int z_ind) const;
	void init_randomwalk_cdf(Lua* lua);
	void re_emit(LorentzHelper *lh, const int z_ind) const;
	void window(LorentzHelper *lh, const int z_ind);
	void window(Particle *p, const int z_ind);

	// solve for temperature and Ye (if steady_state)
	double equilibrium_damping;
	int    equilibrium_itmax;
	double equilibrium_tolerance;
	void   solve_eq_zone_values();
	void   normalize_radiative_quantities();
	double brent_method(const int zone_index, double (*eq_function)(double, void*), const double min, const double max);

	// update temperature and Ye (if !steady_state)
	void update_zone_quantities();

	// stored minimum and maximum values to assure safety
	int max_particles;

	// simulation parameters
	double step_size;
	int    do_annihilation;
	int    radiative_eq;
	int    rank0;
	int    exponential_decay;

	// random walk parameters
	CDFArray randomwalk_diffusion_time;
	LocateArray randomwalk_xaxis;
	double randomwalk_sphere_size;
	double randomwalk_min_optical_depth;
	double randomwalk_max_x;
	int randomwalk_sumN;
	int randomwalk_n_isotropic;

	// output parameters
	int write_zones_every;
	int write_rays_every;
	int write_spectra_every;

	// global radiation quantities
	double particle_total_energy;
	double particle_rouletted_energy;
	double particle_core_abs_energy;
	double particle_escape_energy;

	// check parameters
	void check_parameters() const;

	// check quality
	int number_of_bins() const;

public:

	Transport();
	virtual ~Transport() {}

	int verbose;

	int    equilibrium_T;
	int    equilibrium_Ye;
	double T_min,  T_max;
	double Ye_min, Ye_max;
	double rho_min, rho_max;

	// arrays of species
	std::vector<Species*> species_list;

	// pointer to grid
	Grid *grid;

	// biasing
	// minimum neutrino packet energy
	double min_packet_energy;
	double max_packet_energy;
	double importance_bias;
	double min_importance;
	int bias_path_length;
	double max_path_length_boost;
	void sample_tau(LorentzHelper *lh);

	// items for core emission
	double r_core;
	int n_emit_core;
	int n_emit_core_per_bin;
	double core_lum_multiplier;
	int core_emit_method;
	CDFArray core_species_luminosity;
	void init_core(const double r_core, const double T_core, const double munue_core);
	void init_core(const double r_core, const std::vector<double>& T_core, const std::vector<double>& mu_core, const std::vector<double>& L_core);

	// items for zone emission
	int do_visc;
	int do_relativity;
	int use_scattering_kernels;
	int n_emit_zones;
	int n_emit_zones_per_bin;
	double visc_specific_heat_rate;

	// how many times do we emit+propagate each timestep?
	int n_subcycles;
	int do_emit_by_bin;

	// global radiation quantities
	std::vector<double> L_core_lab;
	std::vector<double> L_net_lab;
	std::vector<double> L_net_esc;
	std::vector<double> N_net_lab;
	std::vector<double> N_net_esc;
	std::vector<long> n_active;
	std::vector<long> n_escape;


	// random number generator
	mutable ThreadRNG rangen;

	// set things up
	void   init(Lua* lua);

	// in-simulation functions to be used by main
	void step();
	void write(const int it) const;
	int  total_particles() const;
	void write_rays(const int it);
	static std::string filename(const char* filebase, const int iw, const char* suffix);
	static double mean_mass(const double Ye);
	double importance(const double abs_opac, const double scat_opac, const double dx) const;


	// per-zone luminosity functions
	double zone_comoving_visc_heat_rate(const int zone_index) const;
	double zone_comoving_therm_emit_energy (const int zone_index) const;
	double zone_comoving_therm_emit_leptons(const int zone_index) const;
	double zone_comoving_biased_therm_emit_energy(const int z_ind) const;
	double bin_comoving_therm_emit_energy(const int z_ind, const int s, const int g) const;

};

#endif

