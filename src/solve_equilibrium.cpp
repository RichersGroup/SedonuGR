#pragma warning disable 161
#include <omp.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include "transport.h"
#include "Lua.h"
#include "grid_general.h"
#include "species_general.h"
#include "global_options.h"

double temp_eq_function(int z_ind, double T,  transport* sim);
double   Ye_eq_function(int z_ind, double Ye, transport* sim);


//-------------------------------------------------------------
//  Solve for the temperature assuming radiative equilibrium
//-------------------------------------------------------------
void transport::solve_eq_zone_values()
{
	assert(brent_solve_tolerance > 0);

	// remember what zones I'm responsible for
	int start = ( MPI_myID==0 ? 0 : my_zone_end[MPI_myID - 1] );
	int end = my_zone_end[MPI_myID];
	assert(end >= start);
	assert(start >= 0);
	assert(end < grid->z.size());

	// solve radiative equilibrium temperature and Ye (but only in the zones I'm responsible for)
	// don't solve if out of density bounds
#pragma omp parallel for schedule(guided)
	for (int z_ind=start; z_ind<end; z_ind++) if( (grid->z[z_ind].rho >= rho_min) && (grid->z[z_ind].rho <= rho_max) )
	{
		double T_last_iter=NaN, Ye_last_iter=NaN;
		double T_last_step=NaN, Ye_last_step=NaN;
		double T_error=NaN,Ye_error=NaN;
		double dT_step=NaN, dYe_step=NaN;

		// set up the solver
		if(solve_T)
		{
			T_error  = 10*brent_solve_tolerance;
			T_last_step  = grid->z[z_ind].T_gas;
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
				T_last_iter  = grid->z[z_ind].T_gas;
				grid->z[z_ind].T_gas = brent_method(z_ind, temp_eq_function, T_min,  T_max);
				T_error  = fabs( (grid->z[z_ind].T_gas - T_last_iter ) / (T_last_iter ) );
			}
			if(solve_Ye)
			{
				Ye_last_iter = grid->z[z_ind].Ye;
				grid->z[z_ind].Ye    = brent_method(z_ind, Ye_eq_function,   Ye_min, Ye_max);
				Ye_error = fabs( (grid->z[z_ind].Ye    - Ye_last_iter) / (Ye_last_iter) );
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
				dT_step  = grid->z[z_ind].T_gas - T_last_step;
				grid->z[z_ind].T_gas =  T_last_step + (1.0 - damping)*dT_step;
				if(grid->z[z_ind].T_gas > T_max){
					cout << "# WARNING: Changing T_gas in zone " << z_ind << " from " << grid->z[z_ind].T_gas << " to T_max=" << T_max << endl;
					grid->z[z_ind].T_gas = T_max;}
				if(grid->z[z_ind].T_gas < T_min){
					cout << "# WARNING: Changing T_gas in zone " << z_ind << " from " << grid->z[z_ind].T_gas << " to T_min=" << T_min << endl;
					grid->z[z_ind].T_gas = T_min;}
				if(grid->z[z_ind].T_gas != grid->z[z_ind].T_gas){
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
double temp_eq_function(int z_ind, double T, transport* sim)
{
	assert(z_ind >= 0);
	assert(z_ind < (int)sim->grid->z.size());
	assert(T >= 0);

	// total energy absorbed in zone
	double E_absorbed = sim->grid->z[z_ind].e_abs;
	// total energy emitted (to be calculated based on emissivities)
	double E_emitted = 0.;

	// set the zone temperature
	sim->grid->z[z_ind].T_gas = T;

	// include the emission from all species
	for(int i=0; i<sim->species_list.size(); i++)
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
	assert(E_emitted > 0);
	assert(E_absorbed > 0);
	return (E_emitted - E_absorbed);
}


//----------------------------------------------------------------------------
// This is the function that expresses radiative equillibrium
// in a cell (i.e. E_absorbed = E_emitted).  It is used in
// The Brent solver below to determine the temperature such
// that RadEq holds
//----------------------------------------------------------------------------
double Ye_eq_function(int z_ind, double Ye, transport* sim)
{
	assert(z_ind >= 0);
	assert(z_ind < (int)sim->grid->z.size());
	assert(Ye >= 0);
	assert(Ye <= 1);

	// total energy absorbed in zone
	double l_absorbed = sim->grid->z[z_ind].l_abs;
	// total energy emitted (to be calculated)
	double l_emitted = 0.;

	// set the zone temperature
	sim->grid->z[z_ind].Ye = Ye;

	// include the emission from all species
	for(int i=0; i<sim->species_list.size(); i++)
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
	}

	// radiative equillibrium condition: "emission equals absorbtion"
	// return to Brent function to iterate this to zero
	return (l_emitted - l_absorbed);
}






//-----------------------------------------------------------
// Brents method (from Numerical Recipes) to solve 
// non-linear equation for T in rad equillibrium
//-----------------------------------------------------------
// definitions used for temperature solver
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
double transport::brent_method(int z_ind, double (*eq_function)(int,double,transport*), double min, double max)
{
	assert(z_ind >= 0);

	double small = 3.0e-8;
	int iter;

	// Initial guesses
	double a=min;
	double b=max;
	double c=b;
	double d,e,min1,min2;
	double fa=(*eq_function)(z_ind,a,this);
	double fb=(*eq_function)(z_ind,b,this);
	double fc,p,q,r,s,tol1,xm;

	//if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0))
	//  printf("Root must be bracketed in zbrent");
	fc=fb;
	for (iter=1;iter<=brent_itmax;iter++) {
		if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
			c=a;
			fc=fa;
			e=d=b-a;
		}
		if (fabs(fc) < fabs(fb)) {
			a=b;
			b=c;
			c=a;
			fa=fb;
			fb=fc;
			fc=fa;
		}
		tol1=2.0*small*fabs(b)+0.5*brent_solve_tolerance;
		xm=0.5*(c-b);
		if (fabs(xm) <= tol1 || fb == 0.0) return b;
		if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
			s=fb/fa;
			if (a == c) {
				p=2.0*xm*s;
				q=1.0-s;
			} else {
				q=fa/fc;
				r=fb/fc;
				p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
				q=(q-1.0)*(r-1.0)*(s-1.0);
			}
			if (p > 0.0) q = -q;
			p=fabs(p);
			min1=3.0*xm*q-fabs(tol1*q);
			min2=fabs(e*q);
			if (2.0*p < (min1 < min2 ? min1 : min2)) {
				e=d;
				d=p/q;
			} else {
				d=xm;
				e=d;
			}
		} else {
			d=xm;
			e=d;
		}
		a=b;
		fa=fb;
		if (fabs(d) > tol1)
			b += d;
		else
			b += SIGN(tol1,xm);
		fb=(*eq_function)(z_ind,b,this);
	}
	printf("Maximum number of iterations exceeded in zbrent\n");
	return 0.0;
}
