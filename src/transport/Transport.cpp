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

#include <omp.h>
#include <mpi.h>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <string>
#include "physical_constants.h"
#include "Transport.h"
#include "Grid.h"
#include "Grid0DIsotropic.h"
#include "Grid1DSphere.h"
#include "Grid2DSphere.h"
#include "Grid3DCart.h"
#include "GridGR1D.h"
#include "Species.h"
#include "Neutrino_NuLib.h"
#include "Neutrino_GR1D.h"
#include "Neutrino_Nagakura.h"
#include "Neutrino_grey.h"
#include "nulib_interface.h"
#include "global_options.h"

using namespace std;
namespace pc = physical_constants;

// constructor
Transport::Transport(){
	verbose = -MAXLIM;
	MPI_nprocs = -MAXLIM;
	MPI_myID = -MAXLIM;
	T_min = NaN;
	T_max = NaN;
	Ye_min = NaN;
	Ye_max = NaN;
	rho_min = NaN;
	rho_max = NaN;
	min_step_size = NaN;
	max_step_size = NaN;
	do_randomwalk = -MAXLIM;
	min_packet_weight = NaN;
	do_annihilation = -MAXLIM;
	grid = NULL;
	r_core = NaN;
	n_emit_core_per_bin = -MAXLIM;
	n_emit_zones_per_bin = -MAXLIM;
	n_subcycles = -MAXLIM;
	write_zones_every = -MAXLIM;
	particle_core_abs_energy = NaN;
	particle_rouletted_energy = NaN;
	particle_escape_energy = NaN;
	randomwalk_min_optical_depth = NaN;
	randomwalk_max_x = NaN;
	randomwalk_sumN = -MAXLIM;
}


//----------------------------------------------------------------------------
// Initialize the transport module
// Includes setting up the grid, particles,
// and MPI work distribution
//----------------------------------------------------------------------------
void Transport::init(Lua* lua)
{ 
	// get mpi rank
	MPI_Comm_size( MPI_COMM_WORLD, &MPI_nprocs );
	MPI_Comm_rank( MPI_COMM_WORLD, &MPI_myID  );
	if(MPI_myID==0){
		cout << "# NDIMS=" << NDIMS << endl;
		cout << "# DEBUG=" << DEBUG << endl;
		cout << "# DO_GR=" << DO_GR << endl;
		cout << "# Initializing transport..." << endl << flush;
		cout << "#   Using " << MPI_nprocs << " MPI ranks" << endl << flush;
#ifdef _OPENMP
#pragma omp parallel
#pragma omp single
		cout << "#   Using " << omp_get_num_threads()  << " threads on each MPI rank." << endl << flush;
#endif
	}

	// figure out what emission models we're using
	n_subcycles = lua->scalar<int>("n_subcycles");
	PRINT_ASSERT(n_subcycles,>=,1);
	n_emit_zones_per_bin = lua->scalar<int>("n_emit_therm_per_bin");
	n_emit_core_per_bin  = lua->scalar<int>("n_emit_core_per_bin");

	// read simulation parameters
	verbose      = MPI_myID==0 ? lua->scalar<int>("verbose") : 0;
	do_annihilation = lua->scalar<int>("do_annihilation");
	min_step_size     = lua->scalar<double>("min_step_size");
	max_step_size     = lua->scalar<double>("max_step_size");
	do_randomwalk = lua->scalar<int>("do_randomwalk");
	if(do_randomwalk){
		randomwalk_min_optical_depth = lua->scalar<double>("randomwalk_min_optical_depth");
		init_randomwalk_cdf(lua);
	}
	min_packet_weight = lua->scalar<double>("min_packet_weight");

	// output parameters
	write_zones_every   = lua->scalar<double>("write_zones_every");


	//===================================//
    // set neutrino's min and max values //
	//===================================//
    T_min   =  0.0;
    T_max   =  numeric_limits<double>::infinity();
    Ye_min  = 0;
    Ye_max  = 1;
    rho_min = -numeric_limits<double>::infinity();
    rho_max =  numeric_limits<double>::infinity();

	// set up nulib table
	string neutrino_type = lua->scalar<string>("neutrino_type");
	if(neutrino_type=="NuLib"){ // get everything from NuLib
		// read the fortran module into memory
	        if(verbose) cout << "# Initializing NuLib..." << endl << flush;
		string nulib_table = lua->scalar<string>("nulib_table");
		nulib_init(nulib_table);

		// eos
		string eos_filename = lua->scalar<string>("nulib_eos");
		nulib_eos_read_table((char*)eos_filename.c_str());

		// set neutrino's min and max values
		T_min  =  nulib_get_Tmin();
		T_max  =  nulib_get_Tmax();
		Ye_min =  nulib_get_Yemin();
		Ye_max =  nulib_get_Yemax();
		rho_min = nulib_get_rhomin();
		rho_max = nulib_get_rhomax();
	}

	//==================//
	// SET UP TRANSPORT //
	//==================//
	int num_nut_species = 0;
	if(neutrino_type=="NuLib"){ // get everything from NuLib
		if(verbose) cout << "#   Using NuLib neutrino opacities" << endl;
		num_nut_species = nulib_get_nspecies();
	}
	else if(neutrino_type=="grey"){
		if(verbose) cout << "#   Using grey opacity (0 chemical potential blackbody)" << endl;
		num_nut_species = 1;
	}
	else if(neutrino_type=="Nagakura"){
		if(verbose) cout << "#   Using Nagakura neutrino opacities, assumed to match the grid" << endl;
		num_nut_species = 3;
	}
	else if(neutrino_type=="GR1D"){
		if(verbose) cout << "#   Using GR1D neutrino opacities, assumed to match the grid" << endl;
		num_nut_species = 3;
	}
	else{
		cout << "ERROR: invalid neutrino type" << endl;
		assert(false);
	}

	// create the species arrays
	if(verbose) cout << "# Setting up misc. transport tools..." << endl << flush;
	for(int i=0; i<num_nut_species; i++){
		Species* neutrinos_tmp;
		if     (neutrino_type == "NuLib")    neutrinos_tmp = new Neutrino_NuLib;
		else if(neutrino_type == "grey" )    neutrinos_tmp = new Neutrino_grey;
		else if(neutrino_type == "Nagakura") neutrinos_tmp = new Neutrino_Nagakura;
		else if(neutrino_type == "GR1D")     neutrinos_tmp = new Neutrino_GR1D;
		else exit(1);
		neutrinos_tmp->ID = i;
		neutrinos_tmp->num_species = num_nut_species;
		species_list.push_back(neutrinos_tmp);
	}


	//=================//
	// SET UP THE GRID //
	//=================//
	// read the grid type
	string grid_type = lua->scalar<string>("grid_type");

	// create a grid of the appropriate type
	if     (grid_type == "Grid0DIsotropic"    ) grid = new Grid0DIsotropic;
	else if(grid_type == "Grid1DSphere"       ) grid = new Grid1DSphere;
	else if(grid_type == "Grid2DSphere"       ) grid = new Grid2DSphere;
	else if(grid_type == "Grid3DCart"         ) grid = new Grid3DCart;
	else if(grid_type == "GridGR1D"           ) {if(verbose) cout << "#   Using GridGR1D" << endl;} // already set up in GR1Dinterface.cpp
	else{
		if(verbose) std::cout << "# ERROR: the requested grid type is not implemented." << std::endl;
		exit(3);}
	grid->init(lua, this);

	//===============//
	// GENERAL SETUP //
	//===============//
	// figure out which zones are in this processors work load
	// a processor will do work in range [start,end)
	my_zone_end.resize(MPI_nprocs);
	for(int proc=0; proc<MPI_nprocs; proc++){
		// how much work does this processor do?
		int my_job = (int)(grid->rho.size()/(1.0*MPI_nprocs));
		if(my_job < 1) my_job = 1;

		// where does this processor start and stop its work? (only the end needs to be stored)
		int my_zone_start = proc*my_job;
		my_zone_end[proc] = my_zone_start + my_job;

		// make sure last guy finishes it all
		if(proc == MPI_nprocs-1) my_zone_end[proc] = grid->rho.size();

		// make sure nobody goes overboard
		if(my_zone_end[proc] >= grid->rho.size()) my_zone_end[proc] = grid->rho.size();
	}

	// setup and seed random number generator(s)
	rangen.init();


	//==========================//
	// INITIALIZE THE NEUTRINOS //
	//==========================//
	for(size_t i=0; i<species_list.size(); i++) species_list[i]->init(lua);

	// complain if we're not simulating anything
	n_active.resize(species_list.size(),0);
	n_escape.resize(species_list.size(),0);
	if(species_list.size() == 0)
	{
		if(MPI_myID==0) cout << "ERROR: you must simulate at least one species of particle." << endl;
		exit(7);
	}

	PRINT_ASSERT(T_min,>=,0);
	PRINT_ASSERT(T_max,>,T_min);
	PRINT_ASSERT(Ye_min,>=,0);
	PRINT_ASSERT(Ye_max,>,Ye_min);
	PRINT_ASSERT(Ye_max,<=,1.0);

	//=================//
	// SET UP THE CORE //
	//=================//
	if(verbose) cout << "# Initializing the core..." << flush;
	r_core = lua->scalar<double>("r_core");   // cm
	if(n_emit_core_per_bin>0){
		vector<double> core_lum_multiplier = lua->vector<double>("core_lum_multiplier");
		vector<double> T_core = lua->vector<double>("T_core");
		vector<double> mu_core = lua->vector<double>("core_chem_pot");
		PRINT_ASSERT(core_lum_multiplier.size(),==,species_list.size());
		PRINT_ASSERT(T_core.size(),==,species_list.size());
		PRINT_ASSERT(mu_core.size(),==,species_list.size());
		for(size_t s=0; s<species_list.size(); s++){
			species_list[s]->core_lum_multiplier = core_lum_multiplier[s];
			species_list[s]->T_core = T_core[s] / pc::k_MeV;    // K;
			species_list[s]->mu_core = mu_core[s] * pc::MeV_to_ergs; // erg;
		}
	}
	if(verbose) cout << "finished." << endl << flush;

	// check the parameters
	if(verbose) cout << "# Checking parameters..." << flush;
	check_parameters();
	if(verbose) cout << "finished." << endl << flush;

	// explicitly set global radiation quantities to 0
	N_core_emit.resize(species_list.size());
	L_net_esc.resize(species_list.size());
	N_net_emit.resize(species_list.size());
	N_net_esc.resize(species_list.size());
	for(size_t s=0; s<species_list.size(); s++){
		N_core_emit[s] = 0;
		L_net_esc[s] = 0;
		N_net_emit[s] = 0;
		N_net_esc[s] = 0;
	}
}

void Transport::check_parameters() const{
	if(verbose && do_randomwalk)
		cout << "WARNING: Assumptions in random walk approximation are incompatible with inelastic scattering." << endl;
}

//------------------------------------------------------------
// take a transport time step 
//------------------------------------------------------------
void Transport::step()
{
	// reset radiation quantities
	if(verbose) cout << "# Clearing radiation..." << endl;
	reset_radiation();

	// emit, propagate, and normalize. steady_state means no propagation time limit.
	for(int i=0; i<n_subcycles; i++){
	  if(verbose) cout << "# === Subcycle " << i+1 << "/" << n_subcycles << " ===" << endl;
		emit_particles();
		propagate_particles();
	}
	if(MPI_nprocs>1) sum_to_proc0();      // so each processor has necessary info to solve its zones
	normalize_radiative_quantities();

	// calculate annihilation rates
	if(do_annihilation) calculate_annihilation();
}


//---------------------------------
// write all the necessary output
//---------------------------------
void Transport::write(const int it) const{
	if(MPI_myID==0){
		// write zone state when appropriate
		if(write_zones_every>0) if(it%write_zones_every==0 && write_zones_every>0 && it>0){
			if(verbose) cout << "# writing zone file " << it << endl;
			grid->write_zones(it);
		}
	}
}

//----------------------------
// reset radiation quantities
//------------------------------
void Transport::reset_radiation(){
	// clear global radiation quantities
	for(size_t i=0; i<species_list.size(); i++){
		grid->distribution[i]->wipe();
		grid->spectrum[i].wipe();
		n_active[i] = 0;
		n_escape[i] = 0;
		N_core_emit[i] = 0;
		L_net_esc[i] = 0;
		N_net_emit[i] = 0;
		N_net_esc[i] = 0;
	}

	grid->l_abs.wipe();
	grid->l_emit.wipe();
	grid->fourforce_abs.wipe();
	grid->fourforce_emit.wipe();
	grid->fourforce_annihil.wipe();

	particle_core_abs_energy = 0;
	particle_rouletted_energy = 0;
	particle_escape_energy = 0;

	if(verbose) cout << "# Setting zone transport quantities" << endl << flush;
	for(size_t s=0; s<species_list.size(); s++){
		#pragma omp parallel for
		for(size_t z_ind=0;z_ind<grid->rho.size();z_ind++)
			species_list[s]->set_eas(z_ind,grid);
	}
}

//-----------------------------
// calculate annihilation rates
//-----------------------------
void Transport::calculate_annihilation(){
	grid->fourforce_annihil.mpi_gather(my_zone_end);
	if(verbose) cout << "# Calculating annihilation rates..." << flush;

	// remember what zones I'm responsible for
	if(MPI_myID != 0) return; // fast enough not to have to parallelize
	int start = 0;//( MPI_myID==0 ? 0 : my_zone_end[MPI_myID - 1] );
	int end = grid->rho.size();//my_zone_end[MPI_myID];
	PRINT_ASSERT(end,>=,start);
	PRINT_ASSERT(start,>=,0);
	PRINT_ASSERT(end,<=,(int)grid->rho.size());

	double H_nunu_tet = 0;

    #pragma omp parallel for reduction(+:H_nunu_tet)
	for(int z_ind=start; z_ind<end; z_ind++){

		// get the directional indices
		size_t dir_ind[NDIMS];
		grid->rho.indices(z_ind,dir_ind);

		// get the kernels
		vector< vector< vector< vector<double> > > > phi; // [s][order][gin][gout]
		phi.resize(species_list.size());
		for(size_t s=0; s<species_list.size(); s++)
			species_list[s]->get_annihil_kernels(grid->rho[z_ind], grid->T[z_ind], grid->Ye[z_ind], grid->nu_grid_axis, phi[s]);

		// get the list of species
		vector<Tuple<size_t,2> > pairs;
		switch(species_list.size()){
		case 2:
			pairs.resize(1);
			pairs[0][0]=0; pairs[0][1]=1;
			break;
		case 3:
			pairs.resize(2);
			pairs[0][0]=0; pairs[0][1]=1;
			pairs[1][0]=2; pairs[1][1]=2;
			break;
		case 4:
			pairs.resize(2);
			pairs[0][0]=0; pairs[0][1]=1;
			pairs[1][0]=2; pairs[1][1]=3;
			break;
		case 6:
			pairs.resize(3);
			pairs[0][0]=0; pairs[0][1]=1;
			pairs[1][0]=2; pairs[1][1]=3;
			pairs[2][0]=4; pairs[2][1]=5;
			break;
		default:
			assert(0); // these should be the only options
		}

		for(size_t p=0; p<pairs.size(); p++){
			size_t s0=pairs[p][0], s1=pairs[p][1];
			PRINT_ASSERT(species_list[s0]->weight,==,species_list[s1]->weight);

			grid->distribution[s0]->annihilation_rate(dir_ind,
					grid->distribution[s1],
					phi[s0], species_list[s0]->weight,
					grid->fourforce_annihil[z_ind]);
		}
		H_nunu_tet += grid->fourforce_annihil[z_ind][3] * grid->zone_4volume(z_ind);
	}

	// synchronize global quantities between processors
	// If this gets parallelized, should make mpi_sum more efficient by not transmitting entire array
	//grid->fourforce_annihil.mpi_sum();
	//MPI_Allreduce(MPI_IN_PLACE, &H_nunu_tet, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	// write to screen
	if(verbose) {
		cout << "finished." << endl << flush;
		cout << "#   " << H_nunu_tet << " erg/s H_annihil" << endl << flush;
	}
}


//----------------------------------------------------------------------------
// normalize the radiative quantities
//----------------------------------------------------------------------------
void Transport::normalize_radiative_quantities(){
	if(verbose) cout << "# Normalizing Radiative Quantities" << endl;

	// normalize zone quantities
	double inv_multiplier = 1.0/(double)n_subcycles;
    #pragma omp parallel for
	for(size_t z_ind=0;z_ind<grid->rho.size();z_ind++)
	{
	  //double inv_mult_four_vol = inv_multiplier / grid->zone_4volume(z_ind); // Lorentz invariant - same in lab and comoving frames. Assume lab_dt=1.0
	  //PRINT_ASSERT(inv_mult_four_vol,>=,0);

		grid->fourforce_abs[z_ind] *= inv_multiplier;
		grid->fourforce_emit[z_ind] *= inv_multiplier;
		grid->l_abs[z_ind] *= inv_multiplier;
		grid->l_emit[z_ind] *= inv_multiplier;

		// represents *all* species if nux
		size_t dir_ind[NDIMS];
		grid->rho.indices(z_ind,dir_ind);
		for(size_t s=0; s<species_list.size(); s++){
		  grid->distribution[s]->rescale_spatial_point(dir_ind, inv_multiplier);
		}
	}

	// normalize global quantities
	for(size_t s=0; s<species_list.size(); s++){
		grid->spectrum[s].rescale(inv_multiplier);
		N_core_emit[s] *= inv_multiplier;
		L_net_esc[s] *= inv_multiplier;
		N_net_emit[s] *= inv_multiplier;
		N_net_esc[s] *= inv_multiplier;
	}
	particle_core_abs_energy *= inv_multiplier;
	particle_rouletted_energy *= inv_multiplier;
	particle_escape_energy *= inv_multiplier;

	// output useful stuff
	if(verbose){
		int total_active = 0;
		for(size_t i=0; i<species_list.size(); i++) total_active += n_active[i];
		cout << "#   " << total_active << " particles on all ranks" << endl;

		for(size_t i=0; i<species_list.size(); i++){
			double per_esc = (100.0*(double)n_escape[i])/(double)n_active[i];
			cout << "#     --> " << n_escape[i] << "/" << n_active[i] << " " << species_list[i]->name << " escaped. (" << per_esc << "%)" << endl;
		}

		// particle information (all lab-frame)
		cout << "#   " << particle_rouletted_energy+particle_core_abs_energy+particle_escape_energy << " ERG/S TOTAL PARTICLE END ENERGY " << endl; // assume lab_dt=1.0
		cout << "#     --> " << particle_rouletted_energy << " ERG/S TOTAL ROULETTED PARTICLE ENERGY " << endl; // assume lab_dt=1.0
		cout << "#     --> " << particle_core_abs_energy << " ERG/S TOTAL PARTICLE ENERGY ABSORBED BY CORE" << endl; // assume lab_dt=1.0
		cout << "#     --> " << particle_escape_energy << " ERG/S TOTAL ESCAPED PARTICLE ENERGY " << endl; // assume lab_dt=1.0

		if(n_emit_core_per_bin>0){
			cout << "#   { ";
			for(size_t s=0; s<N_core_emit.size(); s++) cout << setw(12) << N_core_emit[s] << "  ";
			cout << "} 1/s N_core" << endl;
		}

		cout << "#   { ";
		for(size_t s=0; s<L_net_esc.size(); s++) cout << setw(12) << L_net_esc[s] << "  ";
		cout << "} erg/s L_esc (lab)" << endl;

		cout << "#   { ";
		for(size_t s=0; s<N_net_esc.size(); s++) cout << setw(12) << L_net_esc[s]/N_net_esc[s]*pc::ergs_to_MeV << "  ";
		cout << "} MeV E_avg_esc (lab)" << endl;

		cout << "#   { ";
		for(size_t s=0; s<N_net_emit.size(); s++) cout << setw(12) << N_net_emit[s] << "  ";
		cout << "} 1/s N_emit (lab)" << endl;

		cout << "#   { ";
		for(size_t s=0; s<N_net_esc.size(); s++) cout << setw(12) << N_net_esc[s] << "  ";
		cout << "} 1/s N_esc (lab)" << endl;
	}
	//calculate blocking factors
	for(size_t s=0; s<species_list.size(); s++){
		if(verbose) cout << "#     Working on fblock for species " << s << endl;
		for(size_t glob_ind=0;glob_ind<grid->scat_opac[s].size();glob_ind++){
			size_t dir_ind[NDIMS+1];
			grid->scat_opac[s].indices(glob_ind,dir_ind);
			grid->fblock[s][glob_ind]=0.5*(grid->fblock[s][glob_ind]+grid->distribution[s]->return_blocking(dir_ind, species_list[s]->weight));
		}
	}
}


void Transport::sum_to_proc0()
{
	if(verbose) cout << "# Reducing Radiation" << endl;

	// scalars
	if(verbose) cout << "#   Summing scalars" << endl;
	if(MPI_myID==0){
		MPI_Reduce(MPI_IN_PLACE, &particle_rouletted_energy, 1,                  MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(MPI_IN_PLACE, &particle_core_abs_energy,  1,                  MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(MPI_IN_PLACE, &particle_escape_energy,    1,                  MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(MPI_IN_PLACE, &L_net_esc.front(),         L_net_esc.size(),   MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(MPI_IN_PLACE, &N_net_esc.front(),         N_net_esc.size(),   MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(MPI_IN_PLACE, &N_net_emit.front(),        N_net_emit.size(),  MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(MPI_IN_PLACE, &N_core_emit.front(),       N_core_emit.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(MPI_IN_PLACE, &n_escape.front(),          n_escape.size(),    MPI_LONG,   MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(MPI_IN_PLACE, &n_active.front(),          n_active.size(),    MPI_LONG,   MPI_SUM, 0, MPI_COMM_WORLD);
	}
	else{
		MPI_Reduce(&particle_rouletted_energy, NULL,                  1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(&particle_core_abs_energy,  NULL,                  1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(&particle_escape_energy,    NULL,                  1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(&L_net_esc.front(),         NULL,   L_net_esc.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(&N_net_esc.front(),         NULL,   N_net_esc.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(&N_net_emit.front(),        NULL,  N_net_emit.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(&N_core_emit.front(),       NULL, N_core_emit.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(&n_escape.front(),          NULL,    n_escape.size(), MPI_LONG,   MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(&n_active.front(),          NULL,    n_active.size(), MPI_LONG,   MPI_SUM, 0, MPI_COMM_WORLD);
	}
	// volumetric quantities
	if(verbose) cout << "#   Summing interaction rates" << endl;
	grid->fourforce_abs.mpi_sum(); // mpi_sum_scatter(my_zone_end)
	grid->fourforce_emit.mpi_sum();// mpi_sum_scatter(my_zone_end)
	grid->l_abs.mpi_sum();// mpi_sum_scatter(my_zone_end)
	grid->l_emit.mpi_sum();// mpi_sum_scatter(my_zone_end)

	// reduce the spectra and distribution functions
	if(verbose) cout << "#   Summing distribution function & spectra to needed proc and 0" << endl;
	for(size_t i=0; i<species_list.size(); i++){
		if(verbose) cout << "#     Working on species " << i << endl;
		grid->spectrum[i].mpi_sum();
		grid->distribution[i]->mpi_sum();
	}
}

string Transport::filename(const char* filebase, const int iw, const char* suffix){
	string number_string;
	stringstream iwstream;
	iwstream << iw;
	if     (iw < 10)    number_string = "0000" + iwstream.str();
	else if(iw < 100)   number_string = "000"  + iwstream.str();
	else if(iw < 1000)  number_string = "00"   + iwstream.str();
	else if(iw < 10000) number_string = "0"    + iwstream.str();
	else                number_string =          iwstream.str();

	string filename = string(filebase) + "_" + number_string + suffix;
	return filename;
}

double Transport::mean_mass(const double Ye){
	return 1.0 / (Ye/pc::m_p + (1.0-Ye)/pc::m_n);
}

//-----------------------------------------------------------------
// Calculate the fermi-dirac blackbody function (#/s/cm^2/ster/(Hz^3/3))
//-----------------------------------------------------------------
double Transport::number_blackbody(const double T /*K*/, const double chem_pot /*erg*/, const double nu /*Hz*/){
	PRINT_ASSERT(T,>=,0);
	PRINT_ASSERT(chem_pot,==,chem_pot);
	PRINT_ASSERT(nu,>,0);
	if(T==0) return 0;
	double zeta = (pc::h*nu - chem_pot)/pc::k/T;
	double bb = pc::inv_c*pc::inv_c / (exp(zeta) + 1.0);
	PRINT_ASSERT(bb,>=,0);
	return bb;
}


//-----------------------------------------------------
// set cdf to blackbody distribution
// units of emis.N: #/s/cm^2/ster
//-----------------------------------------------------
void Transport::set_cdf_to_BB(const double T, const double chempot, CDFArray& emis){
    #pragma omp parallel for ordered
	for(size_t j=0;j<grid->nu_grid_axis.size();j++)
	{
		double nu  = grid->nu_grid_axis.mid[j];
		double dnu = grid->nu_grid_axis.delta(j);
        #pragma omp ordered
		emis.set_value(j, Transport::number_blackbody(T,chempot,nu)*dnu);
	}
	emis.normalize();
}

void Transport::update_eh_background(EinsteinHelper* eh) const{ // things that depend only on particle position
	// zone index
	eh->z_ind = grid->zone_index(eh->xup);

	// boundary conditions
	if(r_core>0 && radius(eh->xup)<r_core){
	  eh->fate = absorbed;
	  return;
	}
	else if(eh->z_ind<0){
		grid->symmetry_boundaries(eh);
		eh->z_ind = grid->zone_index(eh->xup);
		if(eh->z_ind < 0){
			eh->fate = escaped;
			return;
		}
	}
	grid->grid_coordinates(eh->xup,eh->grid_coords);

	// spatial indices
	grid->rho.indices(eh->z_ind, eh->dir_ind);
	for(size_t i=0; i<NDIMS; i++)	PRINT_ASSERT(eh->dir_ind[i],<,grid->rho.axes[i].size());
	grid->rho.set_InterpolationCube(&(eh->icube_vol),eh->grid_coords,eh->dir_ind);
	eh->icube_vol.set_slope_weights(eh->grid_coords);

	// metric and its derivatives
	if(DO_GR){
		grid->interpolate_metric(eh);
		if(eh->g.gtt >= 0){
			eh->z_ind = -1;
			eh->fate = absorbed;
			return;
		}
	}
	eh->zone_fourvolume = grid->zone_coord_volume(eh->z_ind) * (DO_GR ? eh->g.alpha*sqrt(eh->g.gammalow.det()) : 1.); // ccm*s, assumes dt=1s.
 
	// four-velocity
	eh->v = grid->interpolate_fluid_velocity(*eh);
	eh->set_fourvel();

	// set tetrad
	eh->set_tetrad_basis(grid->tetrad_rotation);
}

// make sure kup is consistent with the new background
// interpolate reaction rates
void Transport::update_eh_k_opac(EinsteinHelper* eh) const{
	if(eh->kup[3] <= 0){
		eh->fate = absorbed;
		eh->z_ind = -1;
		return;
	}

	PRINT_ASSERT(eh->kup,==,eh->kup);
	eh->renormalize_kup();
	eh->grid_coords[NDIMS] = min(eh->nu(), grid->nu_grid_axis.max());
	eh->dir_ind[NDIMS] = min(grid->nu_grid_axis.bin(eh->nu()), (int)grid->nu_grid_axis.size()-1);
	eh->eas_ind = grid->abs_opac[eh->s].direct_index(eh->dir_ind);
	grid->abs_opac[eh->s].set_InterpolationCube(&(eh->icube_spec),eh->grid_coords,eh->dir_ind);
	eh->absopac  =  grid->abs_opac[eh->s].interpolate(eh->icube_spec);
	eh->scatopac = grid->scat_opac[eh->s].interpolate(eh->icube_spec);
	eh->inelastic_scatopac = grid->inelastic_scat_opac[eh->s].interpolate(eh->icube_spec);

	PRINT_ASSERT(eh->absopac,>=,0);
	PRINT_ASSERT(eh->scatopac,>=,0);
	PRINT_ASSERT(eh->inelastic_scatopac,>=,0);
}


// Randomly generate new direction isotropically in comoving frame
void Transport::isotropic_direction(Tuple<double,3>& D, ThreadRNG *rangen){
	double costheta = 2.*rangen->uniform() - 1.;
	double sintheta = sqrt(1. - costheta*costheta);
	double phi = 2.*M_PI*rangen->uniform();
	D[0] = sintheta * cos(phi);
	D[1] = sintheta * sin(phi);
	D[2] = costheta;
	Metric::normalize_Minkowski<3>(D);
}

void Transport::isotropic_kup_tet(Tuple<double,4>& kup_tet, ThreadRNG *rangen){
	PRINT_ASSERT(kup_tet[3],>,0);
	Tuple<double,3> D;
	isotropic_direction(D,rangen);

	kup_tet[0] = kup_tet[3] * D[0];
	kup_tet[1] = kup_tet[3] * D[1];
	kup_tet[2] = kup_tet[3] * D[2];

	PRINT_ASSERT(Metric::dot_Minkowski<4>(kup_tet,kup_tet)/(kup_tet[3]*kup_tet[3]),<,TINY);
}

void Transport::random_core_x(Tuple<double,4>& x) const{
	x[3] = 0;
	Tuple<double,3> x3;
	isotropic_direction(x3,&rangen);

	for(size_t i=0; i<3; i++) x[i] = x3[i] * r_core * (1. + TINY);
	int z_ind = grid->zone_index(x);
	PRINT_ASSERT(z_ind,>=,0);
	PRINT_ASSERT(radius(x),>=,r_core);
}

// given k^x/k_tet^t and u^x and the lab-frame distance to the boundary, return the largest random walk sphere size
double Transport::R_randomwalk(const double kx_kttet, const double ux, const double dlab, const double D) const{
	PRINT_ASSERT(dlab*kx_kttet,>,0); // the displacement and the k vector should be in the same direction
	double b = kx_kttet - ux;
	double a = ux*pc::c * randomwalk_max_x / D;
	double c = -dlab;
	double R = NaN;
	double rad = b*b - 4.*a*c;

	if(rad<0) R = 0; // if no solution, say randomwalk can't be done
	else if(abs(4.*a*c/(b*b)) < sqrt(TINY)){
		R = -c / b;
	}
	else{
		double term1 = -b / (2.*a);
		double term2 = sqrt(rad) / (2.*a);
		R = max(term1 + term2, term1-term2);
	}
	R = max(0.,R);
	return R;
}

//-------------------------------------------------------------
// Sample outgoing neutrino direction and energy
//-------------------------------------------------------------
bool Transport::reject_direction(const double mu, const double delta) const{
	PRINT_ASSERT(mu,>=,-1.);
	PRINT_ASSERT(mu,<=,1.);
	// uniform in direction
	if(delta==0) return false;

	// highly anisotropic - must cut off PDF. Reject if outside bounds
	double min_mu = delta> 1.0 ? delta-2. : -1.0;
	double max_mu = delta<-1.0 ? delta+2. :  1.0;
	if(mu<min_mu || mu>max_mu) return true;

	// If inside bounds, use linear PDF
	double delta_eff = max(min(delta, 1.0), -1.0);
	double mubar = 0.5*(min_mu+max_mu);
	double pdfval = 0.5 + delta_eff * (mu-mubar) / (max_mu-min_mu);
	PRINT_ASSERT(pdfval,<=,1.);
	PRINT_ASSERT(pdfval,>=,0.);
	if(rangen.uniform() > pdfval) return true;
	else return false;
}
