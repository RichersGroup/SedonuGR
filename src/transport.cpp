#include "global_options.h"
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


namespace pc = physical_constants;

// constructor
transport::transport(){
	verbose = -MAXLIM;
	MPI_nprocs = -MAXLIM;
	MPI_myID = -MAXLIM;
	solve_T = -MAXLIM;
	solve_Ye = -MAXLIM;
	damping = NaN;
	brent_itmax = -MAXLIM;
	brent_solve_tolerance = NaN;
	T_min = NaN;
	T_max = NaN;
	Ye_min = NaN;
	Ye_max = NaN;
	rho_min = NaN;
	rho_max = NaN;
	max_particles = -MAXLIM;
	step_size = NaN;
	do_photons = -MAXLIM;
	do_neutrinos = -MAXLIM;
	do_distribution = -MAXLIM;
	iterative = -MAXLIM;
	radiative_eq = -MAXLIM;
	rank0 = -MAXLIM;
	grid = NULL;
	t_now = NaN;
	r_core = NaN;
	n_emit_core = -MAXLIM;
	core_lum_multiplier = NaN;
	do_visc = -MAXLIM;
	n_emit_zones = -MAXLIM;
	visc_specific_heat_rate = NaN;
	L_net_lab = NaN;
	L_esc_lab = NaN;
	reflect_outer = -MAXLIM;
	emissions_per_timestep = -MAXLIM;
	n_initial = -MAXLIM;
	initial_BB_T = NaN;
	initial_BB_munue = NaN;
	ratio_emit_by_bin = NAN;
	write_rays_every = -MAXLIM;
	write_spectra_every = -MAXLIM;
	write_zones_every = -MAXLIM;
}


//----------------------------------------------------------------------------
// Initialize the transport module
// Includes setting up the grid, particles,
// and MPI work distribution
//----------------------------------------------------------------------------
void transport::init(Lua* lua)
{ 
	// get mpi rank
	MPI_Comm_size( MPI_COMM_WORLD, &MPI_nprocs );
	MPI_Comm_rank( MPI_COMM_WORLD, &MPI_myID  );
	rank0 = (MPI_myID==0);
	if(rank0){
		cout << "# Initializing transport..." << endl;
		cout << "#   Using " << MPI_nprocs << " MPI ranks" << endl;
            #ifdef _OPENMP
		    #pragma omp parallel
		    #pragma omp single
		    cout << "#   Using " << omp_get_num_threads()  << " threads on each MPI rank." << endl;
            #endif
	}

	// figure out what emission models we're using
	do_visc      = lua->scalar<int>("do_visc");
	if(do_visc) visc_specific_heat_rate = lua->scalar<double>("visc_specific_heat_rate");
	reflect_outer = lua->scalar<int>("reflect_outer");
	n_emit_core  = lua->scalar<int>("n_emit_core");
	n_emit_zones = lua->scalar<int>("n_emit_therm");
	emissions_per_timestep = lua->scalar<int>("emissions_per_timestep");
	ratio_emit_by_bin = lua->scalar<double>("ratio_emit_by_bin");

	// read simulation parameters
	verbose      = lua->scalar<int>("verbose");
	do_photons   = lua->scalar<int>("do_photons");
	do_neutrinos = lua->scalar<int>("do_neutrinos");
	do_distribution = lua->scalar<int>("do_distribution");
	radiative_eq = lua->scalar<int>("radiative_eq");
	iterative = lua->scalar<int>("iterative");
	solve_T       = lua->scalar<int>("solve_T");
	solve_Ye      = lua->scalar<int>("solve_Ye");
	if(iterative){
		if(solve_T || solve_Ye){
			damping               = lua->scalar<double>("damping");
			brent_itmax           = lua->scalar<int>("brent_itmax");
			brent_solve_tolerance = lua->scalar<double>("brent_tolerance");
		}
	}
	step_size     = lua->scalar<double>("step_size");

	// output parameters
	write_zones_every   = lua->scalar<double>("write_zones_every");
	write_rays_every    = lua->scalar<double>("write_rays_every");
	write_spectra_every = lua->scalar<double>("write_spectra_every");

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
	grid->init(lua);

	// Reserve all the memory we might need right now. Speeds up particle additions.
	max_particles = lua->scalar<int>("max_particles");
	particles.reserve(max_particles + 2*grid->z.size()); // to allow for at most 1 additional particle in each cell

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
	n_active.resize(species_list.size(),0);
	n_escape.resize(species_list.size(),0);
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

	// set up initial particle creation
	n_initial = lua->scalar<int>("n_initial");
	if(n_initial>0){
		initial_BB_T     = lua->scalar<double>("initial_BB_T")/pc::k_MeV; // K
		initial_BB_munue = lua->scalar<double>("initial_BB_munue")*pc::MeV_to_ergs; // erg
	}

	//=================//
	// SET UP THE CORE //
	//=================//
	core_lum_multiplier = lua->scalar<double>("core_lum_multiplier");
	double T_core = lua->scalar<double>("T_core") / pc::k_MeV;    // K
	r_core = lua->scalar<double>("r_core");   // cm
	double munue_core = lua->scalar<double>("core_nue_chem_pot") * pc::MeV_to_ergs; // erg
	int iorder = lua->scalar<int>("cdf_interpolation_order");
	core_species_luminosity.interpolation_order = iorder;
	for(unsigned s=0; s<species_list.size(); s++) species_list[s]->core_emis.interpolation_order = iorder;
	if(n_emit_core>0) init_core(r_core, T_core, munue_core);

	// check the parameters
	check_parameters();
}

//-----------------------------------
// set up core (without reading lua)
//-----------------------------------
void transport::init_core(const double r_core /*cm*/, const double T_core /*K*/, const double munue_core /*erg*/){
	assert(n_emit_core>0);
	assert(r_core>0);
	assert(T_core>0);
	assert(species_list.size()>0);

	// set up core emission spectrum function (erg/s)
	core_species_luminosity.resize(species_list.size());
	for(unsigned s=0; s<species_list.size(); s++){
		double chempot = munue_core * (double)species_list[s]->lepton_number; // erg
		if(n_emit_core > 0) species_list[s]->set_cdf_to_BB(T_core, chempot, species_list[s]->core_emis);
		species_list[s]->core_emis.N *= pc::pi * (4.0*pc::pi*r_core*r_core) * species_list[s]->weight * core_lum_multiplier; // erg/s
		core_species_luminosity.set_value(s, species_list[s]->integrate_core_emis());
	}
	core_species_luminosity.normalize();

	// check emission quality from core
	if(verbose && rank0) cout << "# Q_core = " << Q_core() << endl;
}


// estimate the quality of the core emission (#emitted per bin)
double transport::Q_core() const{
	assert(n_emit_core>0);
	assert(r_core>0);
	assert(species_list.size()>0);

	double min_emis = numeric_limits<double>::infinity();
	double net_emis = 0;
	for(unsigned s=0; s<species_list.size(); s++){
		for(unsigned g=0; g<species_list[s]->core_emis.size(); g++){
			double this_emis = species_list[s]->core_emis.get_value(g);
			if(this_emis < min_emis) min_emis = this_emis;
			net_emis += this_emis;
		}
	}
	return min_emis/net_emis * (double)n_emit_core;
}


// estimate the quality of the zone emission (#emitted per bin)
double transport::Q_zones() const{
	assert(n_emit_zones>0);
	assert(species_list.size()>0);
	assert(grid->z.size()>0);

	double min_emis = numeric_limits<double>::infinity();
	double net_emis = 0;
	for(unsigned s=0; s<species_list.size(); s++){
		for(unsigned z_ind=0; z_ind<grid->z.size(); z_ind++){
			double this_min_emis = species_list[s]->min_bin_emis(z_ind);
			if(this_min_emis < min_emis) min_emis = this_min_emis;
			net_emis += species_list[s]->integrate_zone_emis(z_ind);
		}
	}
	return min_emis/net_emis * (double)n_emit_zones;
}

void transport::check_parameters() const{
	if(n_emit_zones>0 && radiative_eq){
		cout << "ERROR: Emitting particles at beginning of timestep AND re-emitting them is inconsistent." << endl;
		exit(10);
	}
	assert(ratio_emit_by_bin >= 0.0);
	assert(ratio_emit_by_bin <= 1.0);
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
	double emission_time = (lab_dt<0 ? 1.0 : lab_dt);
	if(iterative) assert(particles.empty());

	// distribute initial particles in the simulation area
	// don't initialize them if iterative calculation. They all come from the core.
	if(n_initial>0) initialize_blackbody(initial_BB_T,initial_BB_munue);

	// reset radiation quantities
	reset_radiation();

	// emit, propagate, and normalize. steady_state means no propagation time limit.
	for(int i=0; i<emissions_per_timestep; i++){
	  if(rank0 && verbose) cout << "#   subcycle " << i+1 << "/" << emissions_per_timestep << endl;
		emit_particles(emission_time);
		propagate_particles(lab_dt);
	}
	if(MPI_nprocs>1) reduce_radiation();      // so each processor has necessary info to solve its zones
	normalize_radiative_quantities(emission_time);

	// solve for T_gas and Ye structure
	if(solve_T || solve_Ye){
		if(iterative) solve_eq_zone_values();  // solve T,Ye s.t. E_abs=E_emit and N_abs=N_emit
		else update_zone_quantities();            // update T,Ye based on heat capacity and number of leptons
	}
	if(MPI_nprocs>1) synchronize_gas();       // each processor broadcasts its solved zones to the other processors

	// post-processing
	if(rank0) calculate_timescales();
	if(do_distribution && do_neutrinos) calculate_annihilation();

	// advance time step
	if (!iterative) t_now += lab_dt;
}


//---------------------------------
// write all the necessary output
//---------------------------------
void transport::write(const int it) const{
	if(rank0){
		// write zone state when appropriate
		if(it%write_zones_every==0 && write_zones_every>0){
			if(verbose) cout << "# writing zone file " << it << endl;
			grid->write_zones(it);
		}

		// write ray data when appropriate
		if(it%write_rays_every==0 && write_rays_every>0){
			if(verbose) cout << "# writing ray file " << it << endl;
			grid->write_rays(it);
		}

		// print out spectrum when appropriate
		if(it%write_spectra_every==0 && write_spectra_every>0){
			if(verbose) cout << "# writing spectrum files " << it << endl;
			for(unsigned i=0; i<species_list.size(); i++) species_list[i]->spectrum.print(it,i);
		}
	}
}

//----------------------------
// reset radiation quantities
//------------------------------
void transport::reset_radiation(){
	// clear global radiation quantities
	L_net_lab = 0;
	L_esc_lab = 0;
	for(unsigned i=0; i<species_list.size(); i++){
		species_list[i]->spectrum.wipe();
		n_active[i] = 0;
		n_escape[i] = 0;
	}

	#pragma omp parallel
	{
		// calculate the zone eas variables
        #pragma omp for collapse(2)
		for(unsigned i=0; i<species_list.size(); i++)
			for(unsigned j=0; j<grid->z.size(); j++)
				species_list[i]->set_eas(j);

		// prepare zone quantities for another round of transport
        #pragma omp for
		for(unsigned z_ind=0;z_ind<grid->z.size();z_ind++)
		{
			zone* z = &(grid->z[z_ind]);
			z->e_rad    = 0;
			z->nue_abs  = 0;
			z->anue_abs = 0;
			z->e_abs    = 0;
			z->l_emit   = 0;
			z->e_emit   = 0;
			z->Q_annihil = 0;

			if(do_distribution) for(unsigned s=0; s<species_list.size(); s++) z->distribution[s].wipe();
		}
	} // #pragma omp parallel
}

//-----------------------------
// calculate annihilation rates
//-----------------------------
void transport::calculate_annihilation() const{
	if(rank0 && verbose) cout << "# Calculating annihilation rates...";
	double Q_annihil_net = 0;

    #pragma omp parallel for reduction(+:Q_annihil_net)
	for(unsigned z_ind=0; z_ind<grid->z.size(); z_ind++){

		// based on nulib's prescription for which species each distribution represents
		unsigned s_nue    = 0;
		unsigned s_nubare = 1;
		double Q = neutrinos::annihilation_rate(grid->z[z_ind].distribution[s_nue   ],
												grid->z[z_ind].distribution[s_nubare],
												true);
		if(species_list.size()==3){
			unsigned s_nux = 2;
			Q += 2.0 * neutrinos::annihilation_rate(grid->z[z_ind].distribution[s_nux],
													grid->z[z_ind].distribution[s_nux],
													false);
		}
		else if(species_list.size()==4){
			unsigned s_nux = 2;
			unsigned s_nubarx = 3;
			Q += 2.0 * neutrinos::annihilation_rate(grid->z[z_ind].distribution[s_nux   ],
					  	  	  	  	  	  	  	    grid->z[z_ind].distribution[s_nubarx],
					  	  	  	  	  	  	  	    false);
		}
		else if(species_list.size()==6){
			unsigned s_numu = 2;
			unsigned s_nubarmu = 3;
			unsigned s_nutau = 4;
			unsigned s_nubartau = 5;
			Q += neutrinos::annihilation_rate(grid->z[z_ind].distribution[s_numu   ],
											  grid->z[z_ind].distribution[s_nubarmu],
											  false);
			Q += neutrinos::annihilation_rate(grid->z[z_ind].distribution[s_nutau   ],
											  grid->z[z_ind].distribution[s_nubartau],
											  false);
		}
		else{
			cout << "ERROR: wrong species list size in calculate_annihilation" << endl;
			exit(34);
		}

		// set the value in grid
		assert(Q >= 0);
		grid->z[z_ind].Q_annihil = Q;                      // erg/ccm/s
		Q_annihil_net += Q * grid->zone_lab_volume(z_ind); // erg/s
	}
	if(rank0) cout << "finished." << endl;
	if(rank0 && verbose) cout << "# Net neutrino annihilation rate: " << Q_annihil_net << " erg/s" << endl;
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
		z->t_labs  = z->rho/pc::m_p / (z->nue_abs-z->anue_abs);
		z->t_lemit = z->rho/pc::m_p / z->l_emit;
	}
}

//----------------------------------------------------------------------------
// normalize the radiative quantities
//----------------------------------------------------------------------------
void transport::normalize_radiative_quantities(const double lab_dt){
	if(verbose && rank0) cout << "# Normalizing Radiative Quantities" << endl;
	double net_visc_heating = 0;

	// normalize zone quantities
	double multiplier = (double)emissions_per_timestep;
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
			assert(z->nue_abs == 0.0);
			assert(z->anue_abs == 0.0);
			assert(z->e_emit == 0.0);
			assert(z->l_emit == 0.0);
		}

		z->e_rad    /= multiplier*four_vol*pc::c; // erg*dist --> erg/ccm
		z->e_abs    /= multiplier*four_vol;       // erg      --> erg/ccm/s
		z->e_emit   /= multiplier*four_vol;       // erg      --> erg/ccm/s
		z->nue_abs  /= multiplier*four_vol;       // num      --> num/ccm/s
		z->anue_abs /= multiplier*four_vol;       // num      --> num/ccm/s
		z->l_emit   /= multiplier*four_vol;       // num      --> num/ccm/s

		// erg*dist --> erg/ccm, represents *one* species, not *all* of them
		if(do_distribution) for(unsigned s=0; s<species_list.size(); s++)
			z->distribution[s].rescale(1./(multiplier*four_vol*pc::c * species_list[s]->weight));
	}

	// normalize global quantities
	L_net_lab /= multiplier*lab_dt;
	L_esc_lab /= multiplier*lab_dt;
	for(unsigned s=0; s<species_list.size(); s++) species_list[s]->spectrum.rescale(1./(multiplier*lab_dt * species_list[s]->weight));

	// output useful stuff
	if(rank0 && verbose){
		unsigned long total_active = 0;
		for(unsigned i=0; i<species_list.size(); i++){
			double per_esc = (100.0*(double)n_escape[i])/(double)n_active[i];
			total_active += n_active[i];
			cout << "#   " << n_escape[i] << "/" << n_active[i] << " " << species_list[i]->name << " escaped. (" << per_esc << "%)" << endl;
		}
		cout << "#   " << total_active << " total active particles" << endl;
		if(do_visc) cout << "#   " << net_visc_heating << " erg/s Summed comoving-frame viscous heating: " << endl;
		cout << "#   " << L_net_lab << " erg/s Net lab-frame luminosity from (zones+core+re-emission): " << endl;
		cout << "#   " << L_esc_lab << " erg/s Net lab-frame escape luminosity: " << endl;
	}

	if(verbose && rank0 && n_emit_zones>0) cout << "#   Q_zones = " << Q_zones() << endl;
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
	if(verbose && rank0) cout << "# Reducing Radiation" << endl;
	int my_begin, my_end, size;

	// reduce the spectra
	if(verbose && rank0) cout << "#   spectra" << endl;
	for(unsigned i=0; i<species_list.size(); i++) species_list[i]->spectrum.MPI_average();

	// reduce the global luminosity scalars
	if(verbose && rank0) cout << "#   global scalars" << endl;
	double sendscalar;
	sendscalar = L_esc_lab;
	MPI_Reduce(&sendscalar,&L_esc_lab,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	L_esc_lab /= (double)MPI_nprocs;
	sendscalar = L_net_lab;
	MPI_Reduce(&sendscalar,&L_net_lab,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	L_net_lab /= (double)MPI_nprocs;

	// reduce the numbers of particles active/escaped
	if(verbose && rank0) cout << "#   particle numbers" << endl;
	vector<long> sendint(species_list.size(),0);
	vector<long> receiveint(species_list.size(),0);
	for(unsigned i=0; i<species_list.size(); i++) sendint[i] = n_escape[i];
	MPI_Allreduce(&sendint.front(), &receiveint.front(), n_escape.size(), MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
	for(unsigned i=0; i<species_list.size(); i++) n_escape[i] = receiveint[i];
	for(unsigned i=0; i<species_list.size(); i++) sendint[i] = n_active[i];
	MPI_Allreduce(&sendint.front(), &receiveint.front(), n_active.size(), MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
	for(unsigned i=0; i<species_list.size(); i++) n_active[i] = receiveint[i];

	//-- EACH PROCESSOR GETS THE REDUCTION INFORMATION IT NEEDS
	if(verbose && rank0) cout << "#   fluid" << endl;
	//for(int proc=0; proc<MPI_nprocs; proc++){

		// set the begin and end indices so a process covers range [begin,end)
		my_begin = 0;
		my_end = grid->z.size();

		// set the computation size and create the send/receive vectors
		size = my_end - my_begin;
		vector<double> send(size,0);
		vector<double> receive(size,0);

		// reduce distribution
		if(do_distribution)
			for(unsigned s=0; s<species_list.size(); s++)
				for(int i=my_begin; i<my_end; i++) grid->z[i].distribution[s].MPI_average();

		// reduce e_rad
		for(int i=my_begin; i<my_end; i++) send[i-my_begin] = grid->z[i].e_rad;
		MPI_Allreduce(&send.front(), &receive.front(), size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		for(int i=my_begin; i<my_end; i++) grid->z[i].e_rad = receive[i-my_begin] / (double)MPI_nprocs;

		// reduce e_abs
		for(int i=my_begin; i<my_end; i++) send[i-my_begin] = grid->z[i].e_abs;
		MPI_Allreduce(&send.front(), &receive.front(), size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		for(int i=my_begin; i<my_end; i++) grid->z[i].e_abs = receive[i-my_begin] / (double)MPI_nprocs;

		// reduce nue_abs
		for(int i=my_begin; i<my_end; i++) send[i-my_begin] = grid->z[i].nue_abs;
		MPI_Allreduce(&send.front(), &receive.front(), size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		for(int i=my_begin; i<my_end; i++) grid->z[i].nue_abs = receive[i-my_begin] / (double)MPI_nprocs;

		// reduce anue_abs
		for(int i=my_begin; i<my_end; i++) send[i-my_begin] = grid->z[i].anue_abs;
		MPI_Allreduce(&send.front(), &receive.front(), size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		for(int i=my_begin; i<my_end; i++) grid->z[i].anue_abs = receive[i-my_begin] / (double)MPI_nprocs;

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
	//}
}

// after this, only the processor's chunk of radiation quantities
// is correct, but all gas quantities are correct.
void transport::synchronize_gas()
{
	if(verbose && rank0) cout << "# Synchronizing Gas" << endl;
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
	else                             return grid->z[z_ind].H_vis    * grid->z[z_ind].rho * grid->zone_comoving_volume(z_ind);
}


// update zone quantities based on heat capacity and lepton capacity
void transport::update_zone_quantities(){
	if(verbose && rank0) cout << "# Updating Zone Quantities" << endl;
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
			z->Ye += (z->nue_abs-z->anue_abs - z->l_emit) / Nbary;
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

int transport::number_of_bins() const{
	int number_energy_bins = 0;
	for(unsigned s = 0; s<species_list.size(); s++) number_energy_bins += species_list[s]->number_of_bins();
	return number_energy_bins;
}
