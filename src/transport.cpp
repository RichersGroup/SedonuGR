#include <omp.h>
#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <string>
#include "transport.h"
#include "Lua.h"
#include "grid_general.h"
#include "grid_0D_isotropic.h"
#include "grid_1D_sphere.h"
#include "grid_2D_sphere.h"
#include "grid_3D_cart.h"
#include "species_general.h"
#include "photons.h"
#include "neutrinos.h"
#include "cdf_array.h"
#include "nulib_interface.h"
#include "global_options.h"


namespace pc = physical_constants;

// constructor
transport::transport(){
	verbose = -MAX;
	MPI_nprocs = -MAX;
	MPI_myID = -MAX;
	solve_T = -MAX;
	solve_Ye = -MAX;
	damping = NaN;
	brent_itmax = -MAX;
	brent_solve_tolerance = NaN;
	T_min = NaN;
	T_max = NaN;
	Ye_min = NaN;
	Ye_max = NaN;
	rho_min = NaN;
	rho_max = NaN;
	max_particles = -MAX;
	step_size = NaN;
	do_photons = -MAX;
	do_neutrinos = -MAX;
	steady_state = -MAX;
	radiative_eq = -MAX;
	rank0 = -MAX;
	grid = NULL;
	t_now = NaN;
	r_core = NaN;
	n_emit_core = -MAX;
	core_lum_multiplier = NaN;
	do_visc = -MAX;
	n_emit_zones = -MAX;
	visc_specific_heat_rate = NaN;
	L_net_lab = NaN;
	L_esc_lab = NaN;
	reflect_outer = -MAX;
}


//----------------------------------------------------------------------------
// Initialize the transport module
// Includes setting up the grid, particles,
// and MPI work distribution
//----------------------------------------------------------------------------
void transport::init(Lua* lua)
{ 
	if(verbose && rank0) cout << "# Initializing transport..." << endl;
	// get mpi rank
	MPI_Comm_size( MPI_COMM_WORLD, &MPI_nprocs );
	MPI_Comm_rank( MPI_COMM_WORLD, &MPI_myID  );
	rank0 = (MPI_myID==0);

	// figure out what emission models we're using
	n_emit_core  = lua->scalar<int>("n_emit_core");
	n_emit_zones = lua->scalar<int>("n_emit_therm");
	do_visc      = lua->scalar<int>("do_visc");
	if(do_visc) visc_specific_heat_rate = lua->scalar<double>("visc_specific_heat_rate");
	reflect_outer = lua->scalar<int>("reflect_outer");

	// read simulation parameters
	verbose      = lua->scalar<int>("verbose");
	do_photons   = lua->scalar<int>("do_photons");
	do_neutrinos = lua->scalar<int>("do_neutrinos");
	radiative_eq = lua->scalar<int>("radiative_eq");
	steady_state = lua->scalar<int>("steady_state");
	solve_T       = lua->scalar<int>("solve_T");
	solve_Ye      = lua->scalar<int>("solve_Ye");
	if(steady_state){
		if(solve_T || solve_Ye){
			damping               = lua->scalar<double>("damping");
			brent_itmax           = lua->scalar<int>("brent_itmax");
			brent_solve_tolerance = lua->scalar<double>("brent_tolerance");
		}
	}
	step_size     = lua->scalar<double>("step_size");


	// Reserve all the memory we might need right now. Speeds up particle additions.
	max_particles = lua->scalar<int>("max_particles");
	particles.reserve(max_particles);

	//=================//
	// SET UP THE GRID //
	//=================//
	// read the grid type
	string grid_type = lua->scalar<string>("grid_type");

	// create a grid of the appropriate type
	if     (grid_type == "grid_0D_isotropic") grid = new grid_0D_isotropic;
	else if(grid_type == "grid_1D_sphere"   ) grid = new grid_1D_sphere;
	else if(grid_type == "grid_2D_sphere"   ) grid = new grid_2D_sphere;
	else if(grid_type == "grid_3D_cart"     ) grid = new grid_3D_cart;
	else{
		if(rank0) std::cout << "# ERROR: the requested grid type is not implemented." << std::endl;
		exit(3);}

	// initialize the grid (including reading the model file)
	grid->init(lua);

	// calculate integrated quantities to check
	double total_nonrel_mass = 0.0;
	double total_rest_mass   = 0.0;
	double total_rel_KE      = 0.0;
	double total_nonrel_KE   = 0.0;
	double total_rel_TE      = 0.0;
	double total_nonrel_TE   = 0.0;
    #pragma omp parallel for reduction(+:total_nonrel_mass, total_rest_mass, total_rel_KE, total_nonrel_KE, total_rel_TE, total_nonrel_TE)
	for(unsigned z_ind=0;z_ind<grid->z.size();z_ind++){
		double rest_mass   = grid->z[z_ind].rho * grid->zone_comoving_volume(z_ind);
		assert(rest_mass >= 0);
		double nonrel_mass = grid->z[z_ind].rho * grid->zone_lab_volume(z_ind);
		assert(nonrel_mass >= 0);
		vector<double> r;
		grid->zone_coordinates(z_ind,r);

		//if(grid->z[z_ind].rho > 1.0e8 && r[1] > pc::pi/3.0 && r[1] < pc::pi/2.0){
		total_rest_mass += rest_mass;
		total_nonrel_mass += nonrel_mass;
		total_rel_KE    += (rest_mass>0 ? (lorentz_factor(grid->z[z_ind].v) - 1.0) * rest_mass * pc::c*pc::c : 0);
		total_nonrel_KE += 0.5 * nonrel_mass * grid->zone_speed2(z_ind);
		total_rel_TE    += (rest_mass>0 ? rest_mass   / pc::m_n * pc::k * grid->z[z_ind].T : 0);
		total_nonrel_TE += nonrel_mass / pc::m_n * pc::k * grid->z[z_ind].T;
		//}
	}
	if (rank0){
		cout << "#   mass = " << total_rest_mass << " g (nonrel: " << total_nonrel_mass << " g)" <<endl;
		cout << "#   KE = " << total_rel_KE << " erg (nonrel: " << total_nonrel_KE << " erg)" << endl;
		cout << "#   TE = " << total_rel_TE << " erg (nonrel: " << total_nonrel_TE << " erg)" << endl;
	}

	//===============//
	// GENERAL SETUP //
	//===============//
	// figure out which zones are in this processors work load
	// a processor will do work in range [start,end)
	my_zone_end.resize(MPI_nprocs);
	for(int proc=0; proc<MPI_nprocs; proc++){
		// how much work does this processor do?
		int my_job = (int)(grid->z.size()/(1.0*MPI_nprocs));
		if(my_job < 1) my_job = 1;

		// where does this processor start and stop its work? (only the end needs to be stored)
		int my_zone_start = proc*my_job;
		my_zone_end[proc] = my_zone_start + my_job;

		// make sure last guy finishes it all
		if(proc == MPI_nprocs-1) my_zone_end[proc] = grid->z.size();

		// make sure nobody goes overboard
		if(my_zone_end[proc] >= grid->z.size()) my_zone_end[proc] = grid->z.size();
	}

	// setup and seed random number generator(s)
	rangen.init();

	//==================//
	// SET UP TRANSPORT //
	//==================//
	// start at time 0
	t_now = 0;

	/**********************/
	/**/ if(do_photons) /**/
	/**********************/
	{
		photons* photons_tmp = new photons;
		photons_tmp->init(lua, this);
		species_list.push_back(photons_tmp);
	}
	/************************/
	/**/ if(do_neutrinos) /**/
	/************************/
	{
		double grey_opac = lua->scalar<double>("nut_grey_opacity");
		int num_nut_species = 0;
		if(grey_opac < 0){ // get everything from NuLib
			// read the fortran module into memory
			if(rank0) cout << "# Initializing NuLib...";
			string nulib_table = lua->scalar<string>("nulib_table");
			nulib_init(nulib_table);
			if(rank0) cout << "finished." << endl;
			num_nut_species = nulib_get_nspecies();
		}
		else{
			if(rank0) cout << "#   Using grey opacity for electron anti/neutrinos (0 chemical potential)" << endl;
			num_nut_species = 2;
		}

		// create the species arrays
		for(int i=0; i<num_nut_species; i++){
			neutrinos* neutrinos_tmp = new neutrinos;
			neutrinos_tmp->nulibID = i;
			neutrinos_tmp->num_nut_species = num_nut_species;
			neutrinos_tmp->init(lua, this);
			species_list.push_back(neutrinos_tmp);
		}

	}

	// complain if we're not simulating anything
	if(species_list.size() == 0)
	{
		if(rank0) cout << "ERROR: you must simulate at least one species of particle." << endl;
		exit(7);
	}

	// set global min/max values (make range infinite to catch errors)
	T_min   =  numeric_limits<double>::infinity();
	T_max   = -numeric_limits<double>::infinity();
	Ye_min  =  numeric_limits<double>::infinity();
	Ye_max  = -numeric_limits<double>::infinity();
	rho_min =  numeric_limits<double>::infinity();
	rho_max = -numeric_limits<double>::infinity();
	for(unsigned i=0; i<species_list.size(); i++)
	{
		if(species_list[i]->T_min   < T_min  ) T_min   = species_list[i]->T_min;
		if(species_list[i]->T_max   > T_max  ) T_max   = species_list[i]->T_max;
		if(species_list[i]->Ye_min  < Ye_min ) Ye_min  = species_list[i]->Ye_min;
		if(species_list[i]->Ye_max  > Ye_max ) Ye_max  = species_list[i]->Ye_max;
		if(species_list[i]->rho_min < rho_min) rho_min = species_list[i]->rho_min;
		if(species_list[i]->rho_max > rho_max) rho_max = species_list[i]->rho_max;
	}

	assert(T_min >= 0);
	assert(T_max > T_min);
	assert(Ye_min >= 0);
	assert(Ye_max > Ye_min);
	assert(Ye_max <= 1.0);

	// initialize all the zone eas variables
    #pragma omp parallel for collapse(2)
	for(unsigned i=0; i<species_list.size(); i++)
		for(unsigned j=0; j<grid->z.size(); j++)
			species_list[i]->set_eas(j);

	// scatter initial particles in the simulation area
	// don't initialize them if iterative calculation. They all come from the core.
	// int init_particles = lua->scalar<int>("init_particles");
	// if(!steady_state) initialize_particles(init_particles);


	//=================//
	// SET UP THE CORE //
	//=================//
	// the core temperature is used only in setting its emis vector
	// so it's looked at only in species::myInit()
	r_core = lua->scalar<double>("r_core");
	if(n_emit_core > 0){
		core_lum_multiplier = lua->scalar<double>("core_lum_multiplier");
		core_species_luminosity.resize(species_list.size());
		for(unsigned i=0; i<species_list.size(); i++)
			core_species_luminosity.set_value(i, species_list[i]->integrate_core_emis() * core_lum_multiplier);
		core_species_luminosity.normalize();
	}

	// check the parameters
	check_parameters();
}

void transport::check_parameters() const{
	if(n_emit_zones>0 && radiative_eq){
		cout << "ERROR: Emitting particles at beginning of timestep AND re-emitting them is inconsistent." << endl;
		exit(10);
	}
}

double transport::current_time(){
	return t_now;
}

//------------------------------------------------------------
// take a transport time step 
//------------------------------------------------------------
void transport::step(const double lab_dt)
{
	// assume 1.0 s. of particles were emitted if steady_state
	double emission_time = (steady_state ? 1.0 : lab_dt);

    #pragma omp parallel
	{
		// calculate the zone eas variables
        #pragma omp for collapse(2)
		for(unsigned i=0; i<species_list.size(); i++)
			for(unsigned j=0; j<grid->z.size(); j++)
				species_list[i]->set_eas(j);

		// prepare zone quantities for another round of transport
        #pragma omp single
		{
			L_net_lab = 0;
			L_esc_lab = 0;
		}
        #pragma omp for
		for(unsigned z_ind=0;z_ind<grid->z.size();z_ind++)
		{
			zone* z = &(grid->z[z_ind]);
			z->e_rad    = 0;
			z->l_abs    = 0;
			z->e_abs    = 0;
			z->l_emit   = 0;
			z->e_emit   = 0;
		}
	} // #pragma omp parallel

	// emit, propagate, and normalize. steady_state means no propagation time limit.
	emit_particles(emission_time);
	propagate_particles(emission_time);
	normalize_radiative_quantities(emission_time);

	// solve for T_gas and Ye structure
	if(MPI_nprocs>1) reduce_radiation();      // so each processor has necessary info to solve its zones
	if(solve_T || solve_Ye){
		if(steady_state) solve_eq_zone_values();  // solve T,Ye s.t. E_abs=E_emit and N_abs=N_emit
		else update_zone_quantities();            // update T,Ye based on heat capacity and number of leptons
	}
	if(MPI_nprocs>1) synchronize_gas();       // each processor broadcasts its solved zones to the other processors
	calculate_timescales();

	// advance time step
	if (!steady_state) t_now += lab_dt;
}




//-----------------------------
// calculate various timescales
//-----------------------------
void transport::calculate_timescales() const{
#pragma omp parallel for
	for(unsigned i=0;i<grid->z.size();i++){
		zone *z = &(grid->z[i]);
		double e_gas = z->rho*pc::k*z->T/pc::m_p;  // gas energy density (erg/ccm)
		z->t_eabs  = e_gas / z->e_abs;
		z->t_eemit = e_gas / z->e_emit;
		z->t_labs  = z->rho/pc::m_p / z->l_abs;
		z->t_lemit = z->rho/pc::m_p / z->l_emit;
	}
}

//----------------------------------------------------------------------------
// normalize the radiative quantities
//----------------------------------------------------------------------------
void transport::normalize_radiative_quantities(const double lab_dt){
	double net_visc_heating = 0;

	// normalize zone quantities
    #pragma omp parallel for reduction(+:net_visc_heating)
	for(unsigned z_ind=0;z_ind<grid->z.size();z_ind++)
	{
		zone *z = &(grid->z[z_ind]);
		double four_vol = grid->zone_lab_volume(z_ind) * lab_dt; // Lorentz invariant - same in lab and comoving frames

		// add heat absorbed from viscosity to tally of e_abs
		if(do_visc && grid->good_zone(z_ind)){
			z->e_abs += zone_comoving_visc_heat_rate(z_ind) * comoving_dt(lab_dt,z_ind); // erg, comoving frame
			net_visc_heating += zone_comoving_visc_heat_rate(z_ind);      // erg/s
		}

		if(!grid->good_zone(z_ind)){
			assert(z->e_abs == 0.0);
			assert(z->l_abs == 0.0);
			assert(z->e_emit == 0.0);
			assert(z->l_emit == 0.0);
		}

		z->e_rad    /= four_vol*pc::c; // erg*dist --> erg/ccm
		z->e_abs    /= four_vol;       // erg      --> erg/ccm/s
		z->e_emit   /= four_vol;       // erg      --> erg/ccm/s
		z->l_abs    /= four_vol;       // num      --> num/ccm/s
		z->l_emit   /= four_vol;       // num      --> num/ccm/s
	}

	// normalize global quantities
	L_net_lab /= lab_dt;
	L_esc_lab /= lab_dt;

	// output useful stuff
	if(rank0 && verbose){
		if(do_visc) cout << "#   " << net_visc_heating << " erg/s Summed comoving-frame viscous heating: " << endl;
		cout << "#   " << L_net_lab << " erg/s Net lab-frame luminosity from (zones+core+re-emission): " << endl;
		cout << "#   " << L_esc_lab << " erg/s Net lab-frame escape luminosity: " << endl;
	}
}


//----------------------------------------------------------------------------
// sum up the number of particles in all species
//----------------------------------------------------------------------------
int transport::total_particles() const{
	return particles.size();
}


//----------------------------------------------------------------------------
// randomly sample the nu-integrated emissivities of all
// species to determine the species of a new particle
// emitted from the core
//----------------------------------------------------------------------------
int transport::sample_core_species() const
{
	// randomly sample the species (precomputed CDF)
	double z = rangen.uniform();
	return core_species_luminosity.get_index(z);
}



//----------------------------------------------------------------------------
// randomly sample the nu-integrated emissivities of all
// species to determine the species of a new particle
// emitted from a zone
//----------------------------------------------------------------------------
// note: could store a zone_species_cdf structure in transport,
// but this would use more memory. Here, trading CPU cycles for 
// memory. If we are CPU limited, we could change this
int transport::sample_zone_species(const int zone_index) const
{
	cdf_array species_cdf;
	double integrated_emis;
	species_cdf.resize(species_list.size());

	// set values and normalize
	for(unsigned i=0; i<species_list.size(); i++)
	{
		integrated_emis = species_list[i]->integrate_zone_emis(zone_index);
		species_cdf.set_value(i,integrated_emis);
	}
	species_cdf.normalize();

	// randomly sample the species
	double z = rangen.uniform();
	return species_cdf.get_index(z);
}


//------------------------------------------------------------
// Combine the radiation tallies in all zones
// from all processors using MPI reduce
// after this, radiation quantities on all procs match
//------------------------------------------------------------
void transport::reduce_radiation()
{
	vector<double> send, receive;
	int my_begin, my_end, size;

	//-- EACH PROCESSOR GETS THE REDUCTION INFORMATION IT NEEDS
	for(int proc=0; proc<MPI_nprocs; proc++){

		// set the begin and end indices so a process covers range [begin,end)
		my_begin = 0;
		my_end = grid->z.size();

		// set the computation size and create the send/receive vectors
		size = my_end - my_begin;
		send.resize(size);
		receive.resize(size);

		// reduce e_rad
		for(int i=my_begin; i<my_end; i++) send[i-my_begin] = grid->z[i].e_rad;
		MPI_Allreduce(&send.front(), &receive.front(), size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		for(int i=my_begin; i<my_end; i++) grid->z[i].e_rad = receive[i-my_begin] / (double)MPI_nprocs;

		// reduce e_abs
		for(int i=my_begin; i<my_end; i++) send[i-my_begin] = grid->z[i].e_abs;
		MPI_Allreduce(&send.front(), &receive.front(), size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		for(int i=my_begin; i<my_end; i++) grid->z[i].e_abs = receive[i-my_begin] / (double)MPI_nprocs;

		// reduce l_abs
		for(int i=my_begin; i<my_end; i++) send[i-my_begin] = grid->z[i].l_abs;
		MPI_Allreduce(&send.front(), &receive.front(), size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		for(int i=my_begin; i<my_end; i++) grid->z[i].l_abs = receive[i-my_begin] / (double)MPI_nprocs;

		// reduce e_emit
		for(int i=my_begin; i<my_end; i++) send[i-my_begin] = grid->z[i].e_emit;
		MPI_Allreduce(&send.front(), &receive.front(), size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		for(int i=my_begin; i<my_end; i++) grid->z[i].e_emit = receive[i-my_begin] / (double)MPI_nprocs;

		// reduce l_emit
		for(int i=my_begin; i<my_end; i++) send[i-my_begin] = grid->z[i].l_emit;
		MPI_Allreduce(&send.front(), &receive.front(), size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		for(int i=my_begin; i<my_end; i++) grid->z[i].l_emit = receive[i-my_begin] / (double)MPI_nprocs;

		// format for single reduce
		//my_begin = ( proc==0 ? 0 : my_zone_end[proc-1] );
		//my_end = my_zone_end[proc];
		//for(int i=my_begin; i<my_end; i++) send[i-my_begin] = grid->z[i].l_emit;
		//MPI_Reduce(&send.front(), &receive.front(), size, MPI_DOUBLE, MPI_SUM, proc, MPI_COMM_WORLD);
		//for(int i=my_begin; i<my_end; i++) grid->z[i].l_emit = receive[i-my_begin] / (double)MPI_nprocs;
	}
}

// after this, only the processor's chunk of radiation quantities
// is correct, but all gas quantities are correct.
void transport::synchronize_gas()
{
	vector<double> buffer;
	int my_begin, my_end, size;

	//-- EACH PROCESSOR SENDS THE GRID INFORMATION IT SOLVED
	for(int proc=0; proc<MPI_nprocs; proc++){

		// set the begin and end indices so a process covers range [begin,end)
		my_begin = ( proc==0 ? 0 : my_zone_end[proc-1] );
		my_end = my_zone_end[proc];

		// set the computation size and create the send/receive vectors
		size = my_end - my_begin;
		buffer.resize(size);

		// broadcast T_gas
		if(solve_T){
			if(proc==MPI_myID) for(int i=my_begin; i<my_end; i++) buffer[i-my_begin] = grid->z[i].T;
			MPI_Bcast(&buffer.front(), size, MPI_DOUBLE, proc, MPI_COMM_WORLD);
			if(proc!=MPI_myID) for(int i=my_begin; i<my_end; i++) grid->z[i].T = buffer[i-my_begin];
		}

		// broadcast Ye
		if(solve_Ye){
			if(proc==MPI_myID) for(int i=my_begin; i<my_end; i++) buffer[i-my_begin] = grid->z[i].Ye;
			MPI_Bcast(&buffer.front(), size, MPI_DOUBLE, proc, MPI_COMM_WORLD);
			if(proc!=MPI_myID) for(int i=my_begin; i<my_end; i++) grid->z[i].Ye = buffer[i-my_begin];
		}
	}
}


// rate at which viscosity energizes the fluid (erg/s)
double transport::zone_comoving_visc_heat_rate(const int z_ind) const{
	if(visc_specific_heat_rate >= 0) return visc_specific_heat_rate * grid->z[z_ind].rho * grid->zone_comoving_volume(z_ind);
	else                             return grid->z[z_ind].H_com    * grid->z[z_ind].rho * grid->zone_comoving_volume(z_ind);
}


// Write the spectra to disk (using data since the last write)
void transport::write_spectra(const int it)
{
	for(unsigned i=0; i<species_list.size(); i++){
		if(MPI_nprocs>1) species_list[i]->spectrum.MPI_average();
		species_list[i]->spectrum.print(it,i);
		species_list[i]->spectrum.wipe();
	}
}


// update zone quantities based on heat capacity and lepton capacity
void transport::update_zone_quantities(){
	// remember what zones I'm responsible for
	int start = ( MPI_myID==0 ? 0 : my_zone_end[MPI_myID - 1] );
	int end = my_zone_end[MPI_myID];

	// solve radiative equilibrium temperature and Ye (but only in the zones I'm responsible for)
	// don't solve if out of density bounds
    #pragma omp parallel for schedule(guided)
	for (int i=start; i<end; i++) if( (grid->z[i].rho >= rho_min) && (grid->z[i].rho <= rho_max) )
	{
		zone *z = &(grid->z[i]);

		// adjust the temperature based on the heat capacity (erg/K)
		if(solve_T){
			// assert(z.heat_cap > 0);
			// z.T_gas += (z->e_abs-z->e_emit) / z->heat_cap;
			// assert(z->T_gas >= T_min);
			// assert(z->T_gas <= T_max);
		}

		// adjust the Ye based on the lepton capacity (number of leptons)
		if(solve_Ye){
			double Nbary = grid->zone_rest_mass(i) / mean_mass(i);
			assert(Nbary > 0);
			z->Ye += (z->l_abs-z->l_emit) / Nbary;
			assert(z->Ye >= Ye_min);
			assert(z->Ye <= Ye_max);
		}
	}
}

void transport::open_file(const char* filebase, const int iw, ofstream& outf){
	string number_string;
	if     (iw < 10)    number_string = "0000" + to_string(iw);
	else if(iw < 100)   number_string = "000"  + to_string(iw);
	else if(iw < 1000)  number_string = "00"   + to_string(iw);
	else if(iw < 10000) number_string = "0"    + to_string(iw);
	else                number_string =          to_string(iw);

	string filename = string(filebase) + "_" + number_string + ".dat";
	outf.open(filename.c_str());

}

double transport::mean_mass(const double Ye){
	return 1.0 / (Ye/pc::m_p + (1.0-Ye)/pc::m_n);
}
