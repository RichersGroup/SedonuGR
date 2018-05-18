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
#include "EinsteinHelper.h"

class Species;
class Grid;
enum ParticleEvent {interact, randomwalk, nothing};

class Transport
{

protected:

	// this species' list of particles
	ParticleList particles;

	// MPI stuff
	int MPI_nprocs;
	int MPI_myID;
	void sum_to_proc0();
	std::vector<unsigned> my_zone_end;

	// subroutine for calculating timescales
	void calculate_annihilation();

	// main function to emit particles
	void emit_particles();

	// emit from where?
	void emit_inner_source_by_bin();
	void emit_zones_by_bin();

	// what kind of particle to create?
	void create_surface_particle(ParticleList& output, const double Ep, const unsigned int s, const unsigned int g, const unsigned list_index);
	void create_thermal_particle(ParticleList& output, const int zone_index, const double weight, const unsigned int s, const unsigned int g, const unsigned list_index);

	// propagate the particles
	void propagate_particles();
	void propagate(EinsteinHelper* eh) const;
	virtual void move(EinsteinHelper *eh) const;
	void reset_radiation();
	void which_event(EinsteinHelper* eh, ParticleEvent *event) const;
	void boundary_conditions(EinsteinHelper *eh) const;
	void scatter(EinsteinHelper *eh) const;
	void random_walk(EinsteinHelper *eh) const;
	void init_randomwalk_cdf(Lua* lua);
	void window(EinsteinHelper *eh) const;
	void sample_scattering_final_state(EinsteinHelper* eh, const Tuple<double,4>& kup_tet_old) const;



	// solve for temperature and Ye (if steady_state)
	double equilibrium_damping;
	int    equilibrium_itmax;
	double equilibrium_tolerance;
	void   solve_eq_zone_values();
	void   normalize_radiative_quantities();
	double brent_method(const int zone_index, double (*eq_function)(double, void*), const double min, const double max);

	// simulation parameters
	double min_step_size, max_step_size;
	int    do_annihilation;

	// random walk parameters
	CDFArray randomwalk_diffusion_time;
	Axis randomwalk_xaxis;
	int do_randomwalk;
	double randomwalk_min_optical_depth;
	double randomwalk_max_x;
	int randomwalk_sumN;

	// output parameters
	int write_zones_every;

	// global radiation quantities
	double particle_total_energy;
	double particle_rouletted_energy;
	double particle_core_abs_energy;
	double particle_escape_energy;

	// check parameters
	void check_parameters() const;

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

	// minimum neutrino packet energy
	double min_packet_weight;

	// items for core emission
	double r_core;
	int n_emit_core_per_bin;
	void init_core(const double r_core, const double T_core, const double munue_core);
	void random_core_x(Tuple<double,4>& x) const;

	// items for zone emission
	int do_visc;
	int use_scattering_kernels;
	int n_emit_zones_per_bin;
	double visc_specific_heat_rate;

	// how many times do we emit+propagate each timestep?
	int n_subcycles;

	// global radiation quantities
	std::vector<double> N_core_emit;
	std::vector<double> N_net_emit;
	std::vector<double> N_net_esc;
	std::vector<double> L_net_esc;
	std::vector<long> n_active;
	std::vector<long> n_escape;


	// random number generator
	mutable ThreadRNG rangen;

	// blackbody function (#/cm^2/s/ster/Hz^3)
	static double number_blackbody(const double T, const double chempot, const double nu);
	void set_cdf_to_BB(const double T, const double chempot, CDFArray& emis);
	static void isotropic_kup_tet(const double nu, Tuple<double,4>& kup_tet, ThreadRNG *rangen);
	static void isotropic_direction(Tuple<double,3>& D, ThreadRNG *rangen);
	double R_randomwalk(const double kx_kttet, const double kt_kttet, const double ux, const double dlab, const double D);

	// set things up
	void init(Lua* lua);
	void update_eh_background(EinsteinHelper* eh) const;
	void update_eh_k_opac(EinsteinHelper* eh) const;

	// in-simulation functions to be used by main
	void step();
	void write(const int it) const;
	void write_rays(const int it);
	static std::string filename(const char* filebase, const int iw, const char* suffix);
	static double mean_mass(const double Ye);


	// per-zone luminosity functions
	double zone_comoving_visc_heat_rate(const int zone_index) const;
};

#endif

