#include "global_options.h"
#include <omp.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include "transport.h"
#include "Lua.h"
#include "grid_general.h"
#include "species_general.h"

struct  eq_function_params{int z_ind; transport* sim;};
double    temp_eq_function(double T , void* params);
double  Ye_nue_eq_function(double Ye, void* params);
double Ye_anue_eq_function(double Ye, void* params);
double      Ye_eq_function(double Ye, void* params);


//-------------------------------------------------------------
//  Solve for the temperature assuming radiative equilibrium
//-------------------------------------------------------------
void transport::solve_eq_zone_values()
{
	if(verbose && rank0) cout << "# Solving for equilibrium values" << endl;
	assert(brent_solve_tolerance > 0);

	// remember what zones I'm responsible for
	int start = ( MPI_myID==0 ? 0 : my_zone_end[MPI_myID - 1] );
	int end = my_zone_end[MPI_myID];
	assert(end >= start);
	assert(start >= 0);
	assert(end <= (int)grid->z.size());

	// solve radiative equilibrium temperature and Ye (but only in the zones I'm responsible for)
	// don't solve if out of density bounds
    #pragma omp parallel for schedule(guided)
	for(int z_ind=start; z_ind<end; z_ind++) if( (grid->z[z_ind].rho >= rho_min) && (grid->z[z_ind].rho <= rho_max) )
	{
		double T_last_iter=NaN, Ye_last_iter=NaN;
		double T_last_step=NaN, Ye_last_step=NaN;
		double T_error=NaN,Ye_error=NaN;
		double dT_step=NaN, dYe_step=NaN;

		// set up the solver
		if(solve_T)
		{
			T_error  = 10*brent_solve_tolerance;
			T_last_step  = grid->z[z_ind].T;
		}
		if(solve_Ye)
		{
			Ye_error = 10*brent_solve_tolerance;
			Ye_last_step = grid->z[z_ind].Ye;
		}

		// loop through solving the temperature and Ye until both are within error.
		int iter=0;
		while(iter<=brent_itmax && (T_error>brent_solve_tolerance || Ye_error>brent_solve_tolerance))
		{
			if(solve_T)
			{
				T_last_iter  = grid->z[z_ind].T;
				grid->z[z_ind].T = brent_method(z_ind, temp_eq_function, T_min,  T_max);
				assert(grid->z[z_ind].T >= T_min and grid->z[z_ind].T <= T_max);
				T_error  = fabs( (grid->z[z_ind].T - T_last_iter ) / (T_last_iter ) );
			}
			if(solve_Ye)
			{
				Ye_last_iter = grid->z[z_ind].Ye;
				double nue_abs  = grid->z[z_ind].nue_abs;
				double anue_abs = grid->z[z_ind].anue_abs;
				bool neutrinoless = (nue_abs + anue_abs == 0);
				double weight1 = neutrinoless ? 1.0 : anue_abs / (nue_abs+anue_abs);
				double weight2 = neutrinoless ? 0.0 : nue_abs  / (nue_abs+anue_abs);
				double Ye1 = neutrinoless ?
							 brent_method(z_ind, Ye_eq_function,     Ye_min, Ye_max) :
							 brent_method(z_ind, Ye_nue_eq_function, Ye_min, Ye_max);
				double Ye2 = neutrinoless ? 0.0 :
						     brent_method(z_ind, Ye_anue_eq_function,  Ye_min, Ye_max);
				grid->z[z_ind].Ye = weight1*Ye1 + weight2*Ye2;
				assert(grid->z[z_ind].Ye >= Ye_min and grid->z[z_ind].Ye <= Ye_max);
				Ye_error = fabs( (grid->z[z_ind].Ye - Ye_last_iter) / (Ye_last_iter) );
			}
			iter++;
		}

		// warn if it didn't converge
		if(iter == brent_itmax){
			cout << "# WARNING: outer Brent solver hit maximum iterations. (zone:" << z_ind;
			cout << " processor:" << MPI_myID;
            #ifdef _OPENMP_
			cout << " thread:" << omp_get_thread_num();
            #endif
			cout << ")" << endl;
		}

		// damp the oscillations between steps, ensure that it's within the allowed boundaries
		if(damping>0)
		{
			if(solve_T)
			{
				dT_step  = grid->z[z_ind].T - T_last_step;
				grid->z[z_ind].T =  T_last_step + (1.0 - damping)*dT_step;
				if(grid->z[z_ind].T > T_max){
					cout << "# WARNING: Changing T_gas in zone " << z_ind << " from " << grid->z[z_ind].T << " to T_max=" << T_max << endl;
					grid->z[z_ind].T = T_max;}
				if(grid->z[z_ind].T < T_min){
					cout << "# WARNING: Changing T_gas in zone " << z_ind << " from " << grid->z[z_ind].T << " to T_min=" << T_min << endl;
					grid->z[z_ind].T = T_min;}
				if(grid->z[z_ind].T != grid->z[z_ind].T){
					cout << "# ERROR: T_gas is nan." << endl;
					exit(5);}
			}
			if(solve_Ye)
			{
				dYe_step = grid->z[z_ind].Ye - Ye_last_step;
				grid->z[z_ind].Ye = Ye_last_step + (1.0 - damping)*dYe_step;
				if(grid->z[z_ind].Ye > Ye_max){
					cout << " WARNING: Changing Ye in zone " << z_ind << " from " << grid->z[z_ind].Ye << " to Ye_max=" << Ye_max << endl;
					grid->z[z_ind].Ye = Ye_max;}
				if(grid->z[z_ind].Ye < Ye_min){
					cout << " WARNING: Changing Ye in zone " << z_ind << " from " << grid->z[z_ind].Ye << " to Ye_min=" << Ye_min << endl;
					grid->z[z_ind].Ye = Ye_min;}
				if(grid->z[z_ind].Ye != grid->z[z_ind].Ye){
					cout << "# ERROR: Ye is nan." << endl;
					exit(5);}
			}
		}
	}
}


//----------------------------------------------------------------------------
// This is the function that expresses radiative equillibrium
// in a cell (i.e. E_absorbed = E_emitted).  It is used in
// The Brent solver below to determine the temperature such
// that RadEq holds
//----------------------------------------------------------------------------
double temp_eq_function(double T, void *params)
{
	// read the parameters
	struct eq_function_params *p = (struct eq_function_params *) params;
	const int z_ind = p->z_ind;
	transport* sim = p->sim;

	assert(z_ind >= 0);
	assert(z_ind < (int)sim->grid->z.size());
	assert(T >= 0);
	assert(sim->grid->z[z_ind].e_abs >= 0);

	// total energy absorbed in zone
	double E_absorbed = sim->grid->z[z_ind].e_abs; // + sim->grid->z[z_ind].Q_annihil;
	if(sim->do_visc) E_absorbed += sim->zone_comoving_visc_heat_rate(z_ind) / sim->grid->zone_comoving_volume(z_ind);

	// total energy emitted (to be calculated based on emissivities)
	double E_emitted = 0.;

	// set the zone temperature
	sim->grid->z[z_ind].T = T;

	// include the emission from all species
	for(unsigned i=0; i<sim->species_list.size(); i++)
	{
		// reset the eas variables in this zone
		// TODO OPTIMIZE - only set the emissivity variable
		sim->species_list[i]->set_eas(z_ind);

		// integrate emission over frequency (angle
		// integration gives the 4*PI) to get total
		// radiation energy emitted. Opacities are
		// held constant for this (assumed not to change
		// much from the last time step).
		E_emitted += 4.0*pc::pi * sim->species_list[i]->integrate_zone_emis(z_ind);
	}

	// radiative equillibrium condition: "emission equals absorbtion"
	// return to Brent function to iterate this to zero
	assert(E_emitted >= 0);
	assert(E_absorbed >= 0);
	return (E_emitted - E_absorbed);
}


//----------------------------------------------------------------------------
// This is the function that expresses radiative equillibrium
// in a cell (i.e. E_absorbed = E_emitted).  It is used in
// The Brent solver below to determine the temperature such
// that RadEq holds
//----------------------------------------------------------------------------
double Ye_nue_eq_function(double Ye, void *params)
{
	// read the parameters
	struct eq_function_params *p = (struct eq_function_params *) params;
	const int z_ind = p->z_ind;
	transport* sim = p->sim;

	assert(z_ind >= 0);
	assert(z_ind < (int)sim->grid->z.size());
	assert(Ye >= 0);
	assert(Ye <= 1);

	// total energy absorbed in zone
	double l_absorbed = sim->grid->z[z_ind].nue_abs;
	// total energy emitted (to be calculated)
	double l_emitted = 0.;

	// set the zone temperature
	sim->grid->z[z_ind].Ye = Ye;

	// include the emission from all species
	for(unsigned i=0; i<sim->species_list.size(); i++) if(sim->species_list[i]->lepton_number>0)
	{
		// reset the eas variables in this zone
		// OPTIMIZE - only set the emissivity variable
		sim->species_list[i]->set_eas(z_ind);

		// integrate emissison over frequency (angle
		// integration gives the 4*PI) to get total
		// radiation energy emitted. Opacities are
		// held constant for this (assumed not to change
		// much from the last time step).
		l_emitted += 4.0*pc::pi * sim->species_list[i]->integrate_zone_lepton_emis(z_ind);
		assert(l_emitted>0);
	}

	// radiative equillibrium condition: "emission equals absorbtion"
	// return to Brent function to iterate this to zero
	return (l_emitted - l_absorbed);
}
double Ye_anue_eq_function(double Ye, void *params)
{
	// read the parameters
	struct eq_function_params *p = (struct eq_function_params *) params;
	const int z_ind = p->z_ind;
	transport* sim = p->sim;

	assert(z_ind >= 0);
	assert(z_ind < (int)sim->grid->z.size());
	assert(Ye >= 0);
	assert(Ye <= 1);

	// total energy absorbed in zone
	double l_absorbed = sim->grid->z[z_ind].anue_abs;
	// total energy emitted (to be calculated)
	double l_emitted = 0.;

	// set the zone temperature
	sim->grid->z[z_ind].Ye = Ye;

	// include the emission from all species
	for(unsigned i=0; i<sim->species_list.size(); i++) if(sim->species_list[i]->lepton_number<0)
	{
		// reset the eas variables in this zone
		// OPTIMIZE - only set the emissivity variable
		sim->species_list[i]->set_eas(z_ind);

		// integrate emissison over frequency (angle
		// integration gives the 4*PI) to get total
		// radiation energy emitted. Opacities are
		// held constant for this (assumed not to change
		// much from the last time step).
		// minus sign since integrate_zone_lepton_emis will return a negative number
		l_emitted -= 4.0*pc::pi * sim->species_list[i]->integrate_zone_lepton_emis(z_ind);
		assert(l_emitted>0);
	}

	// radiative equillibrium condition: "emission equals absorbtion"
	// return to Brent function to iterate this to zero
	return (l_emitted - l_absorbed);
}
double Ye_eq_function(double Ye, void *params)
{
	// read the parameters
	struct eq_function_params *p = (struct eq_function_params *) params;
	const int z_ind = p->z_ind;
	transport* sim = p->sim;

	assert(z_ind >= 0);
	assert(z_ind < (int)sim->grid->z.size());
	assert(Ye >= 0);
	assert(Ye <= 1);

	// total energy absorbed in zone
	double l_absorbed = sim->grid->z[z_ind].nue_abs - sim->grid->z[z_ind].anue_abs;
	// total energy emitted (to be calculated)
	double l_emitted = 0.;

	// set the zone temperature
	sim->grid->z[z_ind].Ye = Ye;

	// include the emission from all species
	for(unsigned i=0; i<sim->species_list.size(); i++) if(sim->species_list[i]->lepton_number!=0)
	{
		// reset the eas variables in this zone
		// OPTIMIZE - only set the emissivity variable
		sim->species_list[i]->set_eas(z_ind);

		// integrate emissison over frequency (angle
		// integration gives the 4*PI) to get total
		// radiation energy emitted. Opacities are
		// held constant for this (assumed not to change
		// much from the last time step).
		// minus sign since integrate_zone_lepton_emis will return a negative number
		l_emitted += 4.0*pc::pi * sim->species_list[i]->integrate_zone_lepton_emis(z_ind) * sim->species_list[i]->lepton_number;
	}

	// radiative equillibrium condition: "emission equals absorbtion"
	// return to Brent function to iterate this to zero
	return (l_emitted - l_absorbed);
}





//-----------------------------------------------------------
// Brents method to solve
// non-linear equation for T,Ye in rad equillibrium
//-----------------------------------------------------------
// definitions used for solver
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
double transport::brent_method(int z_ind, double (*eq_function)(double,void*), double min, double max)
{
	assert(z_ind >= 0);

	// check if the root is bracketed
	struct eq_function_params params = {z_ind,this};
	double fa = eq_function(min,&params);
	double fb = eq_function(max,&params);
	if(fa*fb > 0){
		if(fabs(fa)<fabs(fb)) return min;
		else return max;
	}


	// allocate storage for the root solver
	const gsl_root_fsolver_type *T = gsl_root_fsolver_bisection;
	gsl_root_fsolver *s = gsl_root_fsolver_alloc(T);

	// initialize the solver
	gsl_function F;
	F.function = eq_function;
	F.params = &params;
	gsl_root_fsolver_set(s,&F,min,max);

	// do iterative solve
	int status = GSL_CONTINUE;
	int iter = 0;
	while(status == GSL_CONTINUE and iter<brent_itmax){
		status = gsl_root_fsolver_iterate (s);
		double a = gsl_root_fsolver_x_lower(s);
		double b = gsl_root_fsolver_x_upper(s);
		status = gsl_root_test_interval(a,b,0,brent_solve_tolerance);
		iter++;
	}

	// free the memory and return
	double result = gsl_root_fsolver_root(s);
	assert(result >= min);
	assert(result <= max);
	gsl_root_fsolver_free(s);
	return result;
}
