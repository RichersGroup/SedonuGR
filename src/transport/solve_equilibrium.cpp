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

#include <gsl/gsl_roots.h>
#include "Transport.h"
#include "Grid.h"
#include "Species.h"
#include "global_options.h"

using namespace std;
namespace pc = physical_constants;

struct  eq_function_params{int z_ind; Transport* sim;};
double    temp_eq_function(double T , void* params);
double      Ye_eq_function(double Ye, void* params);


//-------------------------------------------------------------
//  Solve for the temperature assuming radiative equilibrium
//-------------------------------------------------------------
void Transport::solve_eq_zone_values()
{
	if(verbose) cout << "# Solving for equilibrium values" << endl;
	PRINT_ASSERT(equilibrium_tolerance,>,0);

	// remember what zones I'm responsible for
	if(MPI_myID != 0) return; // fast enough not to have to parallelize
	int start = 0;//( MPI_myID==0 ? 0 : my_zone_end[MPI_myID - 1] );
	int end = grid->rho.size();//my_zone_end[MPI_myID];
	PRINT_ASSERT(end,>=,start);
	PRINT_ASSERT(start,>=,0);
	PRINT_ASSERT(end,<=,(int)grid->rho.size());

	// solve radiative equilibrium temperature and Ye (but only in the zones I'm responsible for)
	// don't solve if out of density bounds
    #pragma omp parallel for schedule(guided)
	for(int z_ind=start; z_ind<end; z_ind++) if( (grid->rho[z_ind] >= rho_min) && (grid->rho[z_ind] <= rho_max) )
	{
		double T_last_iter=NaN, Ye_last_iter=NaN;
		double T_last_step=NaN, Ye_last_step=NaN;
		double T_error=NaN,Ye_error=NaN;
		double dT_step=NaN, dYe_step=NaN;

		// set up the solver
		if(equilibrium_T)
		{
			T_error  = 10*equilibrium_tolerance;
			T_last_step  = grid->T[z_ind];
		}
		if(equilibrium_Ye)
		{
			Ye_error = 10*equilibrium_tolerance;
			Ye_last_step = grid->Ye[z_ind];
		}

		// loop through solving the temperature and Ye until both are within error.
		int iter=0;
		while(iter<=equilibrium_itmax && (T_error>equilibrium_tolerance || Ye_error>equilibrium_tolerance))
		{
			if(equilibrium_T)
			{
				T_last_iter  = grid->T[z_ind];
				grid->T[z_ind] = brent_method(z_ind, temp_eq_function, T_min,  T_max);
				PRINT_ASSERT(grid->T[z_ind],>=,T_min);
				PRINT_ASSERT(grid->T[z_ind],<=,T_max);
				T_error  = fabs( (grid->T[z_ind] - T_last_iter ) / (T_last_iter ) );
			}
			if(equilibrium_Ye)
			{
				Ye_last_iter = grid->Ye[z_ind];
				grid->Ye[z_ind] = brent_method(z_ind, Ye_eq_function, Ye_min, Ye_max);
				PRINT_ASSERT(grid->Ye[z_ind],>=,Ye_min);
				PRINT_ASSERT(grid->Ye[z_ind],<=,Ye_max);
				Ye_error = fabs( (grid->Ye[z_ind] - Ye_last_iter) / (Ye_last_iter) );
			}
			iter++;
		}

		// warn if it didn't converge
		if(iter == equilibrium_itmax){
			cout << "# WARNING: outer Brent solver hit maximum iterations. (zone:" << z_ind;
			cout << " processor:" << MPI_myID;
            #ifdef _OPENMP_
			cout << " thread:" << omp_get_thread_num();
            #endif
			cout << ")" << endl;
		}

		// damp the oscillations between steps, ensure that it's within the allowed boundaries
		if(equilibrium_damping>0)
		{
			if(equilibrium_T)
			{
				dT_step  = grid->T[z_ind] - T_last_step;
				grid->T[z_ind] =  T_last_step + (1.0 - equilibrium_damping)*dT_step;
				if(grid->T[z_ind] > T_max){
					cout << "# WARNING: Changing T_gas in zone " << z_ind << " from " << grid->T[z_ind] << " to T_max=" << T_max << endl;
					grid->T[z_ind] = T_max;}
				if(grid->T[z_ind] < T_min){
					cout << "# WARNING: Changing T_gas in zone " << z_ind << " from " << grid->T[z_ind] << " to T_min=" << T_min << endl;
					grid->T[z_ind] = T_min;}
				if(grid->T[z_ind] != grid->T[z_ind]){
					cout << "# ERROR: T_gas is nan." << endl;
					exit(5);}
			}
			if(equilibrium_Ye)
			{
				dYe_step = grid->Ye[z_ind] - Ye_last_step;
				grid->Ye[z_ind] = Ye_last_step + (1.0 - equilibrium_damping)*dYe_step;
				if(grid->Ye[z_ind] > Ye_max){
					cout << " WARNING: Changing Ye in zone " << z_ind << " from " << grid->Ye[z_ind] << " to Ye_max=" << Ye_max << endl;
					grid->Ye[z_ind] = Ye_max;}
				if(grid->Ye[z_ind] < Ye_min){
					cout << " WARNING: Changing Ye in zone " << z_ind << " from " << grid->Ye[z_ind] << " to Ye_min=" << Ye_min << endl;
					grid->Ye[z_ind] = Ye_min;}
				if(grid->Ye[z_ind] != grid->Ye[z_ind]){
					cout << "# ERROR: Ye is nan." << endl;
					exit(5);}
			}
		}
	}

	if(equilibrium_Ye) grid->Ye.mpi_gather(my_zone_end);
	if(equilibrium_T) grid->T.mpi_gather(my_zone_end);
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
	Transport* sim = p->sim;

	PRINT_ASSERT(z_ind,>=,0);
	PRINT_ASSERT(z_ind,<,(int)sim->grid->rho.size());
	PRINT_ASSERT(T,>=,0);

	// total energy absorbed in zone
	double E_absorbed = sim->grid->fourforce_abs[z_ind][3]; // + sim->grid->z[z_ind].Q_annihil;
	if(sim->do_visc) E_absorbed += sim->zone_comoving_visc_heat_rate(z_ind) / sim->grid->zone_4volume(z_ind);
	PRINT_ASSERT(E_absorbed,>=,0);

	// total energy emitted (to be calculated based on emissivities)
	double E_emitted = 0.;

	// set the zone temperature
	sim->grid->T[z_ind] = T;

	// set up indices
	unsigned dir_ind[NDIMS+1];
	sim->grid->rho.indices(z_ind,dir_ind);

	// include the emission from all species
	for(unsigned s=0; s<sim->species_list.size(); s++)
	{
		// reset the eas variables in this zone
		sim->species_list[s]->set_eas(z_ind, sim->grid);

		// integrate emission over frequency (angle
		// integration gives the 4*PI) to get total
		// radiation energy emitted. Opacities are
		// held constant for this (assumed not to change
		// much from the last time step).
		// E_emitted += 4.0*pc::pi * sim->species_list[i]->integrate_zone_emis(z_ind);
		for(unsigned g=0; g<sim->grid->nu_grid_axis.size(); g++){
			// get eas global index
			dir_ind[NDIMS] = g;
			unsigned eas_ind = sim->grid->abs_opac[s].direct_index(dir_ind);

			// get this zone's contribution to the emitted leptons
			double tmp = sim->grid->BB[s][eas_ind]/*.interpolate(eh.icube_spec)*/ * sim->grid->abs_opac[s][eas_ind]; // #/s/cm^3/sr/(Hz^3/3)
			tmp *= 4.*pc::pi/*sr*/ * sim->grid->nu_grid_axis.delta3(g)/3.0/*Hz^3/3*/ * pc::h*sim->grid->nu_grid_axis.mid[g]/*erg*/;
			E_emitted += tmp;
		}
	}

	// radiative equillibrium condition: "emission equals absorbtion"
	// return to Brent function to iterate this to zero
	PRINT_ASSERT(E_emitted,>=,0);
	PRINT_ASSERT(E_absorbed,>=,0);
	return (E_emitted - E_absorbed);
}


//----------------------------------------------------------------------------
// This is the function that expresses radiative equillibrium
// in a cell (i.e. E_absorbed = E_emitted).  It is used in
// The Brent solver below to determine the temperature such
// that RadEq holds
//----------------------------------------------------------------------------
double Ye_eq_function(double Ye, void *params)
{
	// read the parameters
	struct eq_function_params *p = (struct eq_function_params *) params;
	const int z_ind = p->z_ind;
	Transport* sim = p->sim;

	PRINT_ASSERT(z_ind,>=,0);
	PRINT_ASSERT(z_ind,<,(int)sim->grid->rho.size());
	PRINT_ASSERT(Ye,>=,0);
	PRINT_ASSERT(Ye,<=,1);

	// total energy absorbed in zone
	double l_absorbed = sim->grid->l_abs[z_ind];
	// total energy emitted (to be calculated)
	double l_emitted = 0.;

	// set the zone temperature
	sim->grid->Ye[z_ind] = Ye;

	// set up indices
	unsigned dir_ind[NDIMS+1];
	sim->grid->rho.indices(z_ind,dir_ind);

	// include the emission from all species
	for(unsigned s=0; s<sim->species_list.size(); s++) if(sim->species_list[s]->lepton_number!=0)
	{
		// reset the eas variables in this zone
		// OPTIMIZE - only set the emissivity variable
		sim->species_list[s]->set_eas(z_ind, sim->grid);

		// integrate emissison over frequency (angle
		// integration gives the 4*PI) to get total
		// radiation energy emitted. Opacities are
		// held constant for this (assumed not to change
		// much from the last time step).
		for(unsigned g=0; g<sim->grid->nu_grid_axis.size(); g++){
			// get eas global index
			dir_ind[NDIMS] = g;
			unsigned eas_ind = sim->grid->abs_opac[s].direct_index(dir_ind);

			// get this zone's contribution to the emitted leptons
			double tmp = sim->grid->BB[s][eas_ind]/*.interpolate(eh.icube_spec)*/ * sim->grid->abs_opac[s][eas_ind]; // #/s/cm^3/sr/(Hz^3/3)
			tmp *= 4.*pc::pi/*sr*/ * sim->grid->nu_grid_axis.delta3(g)/3.0/*Hz^3/3*/ * (float)sim->species_list[s]->lepton_number;
			l_emitted += tmp;
		}
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
double Transport::brent_method(int z_ind, double (*eq_function)(double,void*), double min, double max)
{
	PRINT_ASSERT(z_ind,>=,0);

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
	while(status == GSL_CONTINUE and iter<equilibrium_itmax){
		status = gsl_root_fsolver_iterate (s);
		double a = gsl_root_fsolver_x_lower(s);
		double b = gsl_root_fsolver_x_upper(s);
		status = gsl_root_test_interval(a,b,0,equilibrium_tolerance);
		iter++;
	}

	// free the memory and return
	double result = gsl_root_fsolver_root(s);
	PRINT_ASSERT(result,>=,min);
	PRINT_ASSERT(result,<=,max);
	gsl_root_fsolver_free(s);
	return result;
}
