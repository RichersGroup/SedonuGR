#pragma warning disable 161
#include <omp.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include "transport.h"
#include "physical_constants.h"
#include "Lua.h"
#include "grid_general.h"
#include "species_general.h"

namespace pc = physical_constants;

double temp_eq_function(int zone_index, double T,  transport* sim);
double   Ye_eq_function(int zone_index, double Ye, transport* sim);


//-------------------------------------------------------------
//  Solve for the temperature assuming radiative equilibrium
//-------------------------------------------------------------
// #define brent_solve_tolerance 1.e-2
// #define brent_itmax 100
void transport::solve_eq_zone_values()
{
  // remember what zones I'm responsible for
  int start = ( MPI_myID==0 ? 0 : my_zone_end[MPI_myID - 1] );
  int end = my_zone_end[MPI_myID];

  // solve radiative equilibrium temperature and Ye (but only in the zones I'm responsible for)
  #pragma omp parallel for schedule(guided)
  for (int i=start; i<end; i++)
  {
    double T_last_iter, Ye_last_iter;
    double T_last_step, Ye_last_step;
    double T_error,Ye_error;
    double dT_step, dYe_step;
    int iter=0;

    iter = 0;

    // set up the solver
    if(solve_T)
    {
      T_error  = 10*brent_solve_tolerance;
      T_last_step  = grid->z[i].T_gas;
    }
    if(solve_Ye)
    {
      Ye_error = 10*brent_solve_tolerance;
      Ye_last_step = grid->z[i].Ye;
    }

    // loop through solving the temperature and Ye until both are within error.
    while(iter<=brent_itmax && (T_error>brent_solve_tolerance || Ye_error>brent_solve_tolerance))
    {
      if(solve_T)
      {
	T_last_iter  = grid->z[i].T_gas;
	grid->z[i].T_gas = brent_method(i, temp_eq_function, T_min,  T_max);
	T_error  = fabs( (grid->z[i].T_gas - T_last_iter ) / (T_last_iter ) );
      }
      if(solve_Ye)
      {
	Ye_last_iter = grid->z[i].Ye;
	grid->z[i].Ye    = brent_method(i, Ye_eq_function,   Ye_min, Ye_max);
	Ye_error = fabs( (grid->z[i].Ye    - Ye_last_iter) / (Ye_last_iter) );
      }
      iter++;
    }

    // warn if it didn't converge
    if(iter == brent_itmax) cout << "# WARNING: outer Brent solver hit maximum iterations." << endl;

    // damp the oscillations between steps, ensure that it's within the allowed boundaries
    if(solve_T)
    {
      dT_step  = grid->z[i].T_gas - T_last_step;
      grid->z[i].T_gas =  T_last_step + (1.0 - damping)*dT_step;
      if(grid->z[i].T_gas > T_max){
	cout << "# WARNING: Changing T_gas in zone " << i << " from " << grid->z[i].T_gas << " to T_max=" << T_max << endl;
	grid->z[i].T_gas = T_max;}
      if(grid->z[i].T_gas < T_min){
	cout << "# WARNING: Changing T_gas in zone " << i << " from " << grid->z[i].T_gas << " to T_min=" << T_min << endl;
	grid->z[i].T_gas = T_min;}
    }
    if(solve_Ye)
    {
      dYe_step = grid->z[i].Ye - Ye_last_step;
      grid->z[i].Ye = Ye_last_step + (1.0 - damping)*dYe_step;
      if(grid->z[i].Ye > Ye_max){
	cout << " WARNING: Changing Ye in zone " << i << " from " << grid->z[i].Ye << " to Ye_max=" << Ye_max << endl;
	grid->z[i].Ye = Ye_max;}
      if(grid->z[i].Ye < Ye_min){
	cout << " WARNING: Changing Ye in zone " << i << " from " << grid->z[i].Ye << " to Ye_min=" << Ye_min << endl;
	grid->z[i].Ye = Ye_min;}
    }
  }
}


//----------------------------------------------------------------------------
// This is the function that expresses radiative equillibrium
// in a cell (i.e. E_absorbed = E_emitted).  It is used in
// The Brent solver below to determine the temperature such
// that RadEq holds
//----------------------------------------------------------------------------
double temp_eq_function(int zone_index, double T, transport* sim)
{
  // total energy absorbed in zone
  double E_absorbed = sim->grid->z[zone_index].e_abs;
  // total energy emitted (to be calculated)
  double E_emitted = 0.;

  // set the zone temperature
  sim->grid->z[zone_index].T_gas = T;
  
  // include the emission from all species
  for(int i=0; i<sim->species_list.size(); i++)
  {
    // reset the eas variables in this zone
    // OPTIMIZE - only set the emissivity variable
    sim->species_list[i]->set_eas(zone_index);

    // integrate emisison over frequency (angle
    // integration gives the 4*PI) to get total
    // radiation energy emitted. Opacities are
    // held constant for this (assumed not to change
    // much from the last time step).
    E_emitted += 4.0*pc::pi * sim->species_list[i]->int_zone_emis(zone_index);
  }
  
  // radiative equillibrium condition: "emission equals absorbtion"
  // return to Brent function to iterate this to zero
  return (E_emitted - E_absorbed);
}


//----------------------------------------------------------------------------
// This is the function that expresses radiative equillibrium
// in a cell (i.e. E_absorbed = E_emitted).  It is used in
// The Brent solver below to determine the temperature such
// that RadEq holds
//----------------------------------------------------------------------------
double Ye_eq_function(int zone_index, double Ye, transport* sim)
{
  // total energy absorbed in zone
  double l_absorbed = sim->grid->z[zone_index].l_abs;
  // total energy emitted (to be calculated)
  double l_emitted = 0.;

  // set the zone temperature
  sim->grid->z[zone_index].Ye = Ye;
  
  // include the emission from all species
  for(int i=0; i<sim->species_list.size(); i++)
  {
    // reset the eas variables in this zone
    // OPTIMIZE - only set the emissivity variable
    sim->species_list[i]->set_eas(zone_index);

    // integrate emisison over frequency (angle
    // integration gives the 4*PI) to get total
    // radiation energy emitted. Opacities are
    // held constant for this (assumed not to change
    // much from the last time step).
    l_emitted += 4.0*pc::pi * sim->species_list[i]->int_zone_lepton_emis(zone_index);
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
double transport::brent_method(int zone_index, double (*eq_function)(int,double,transport*), double min, double max)
{
  double small = 3.0e-8;
  int iter;

  // Initial guesses
  double a=min;
  double b=max;
  double c=b;
  double d,e,min1,min2;
  double fa=(*eq_function)(zone_index,a,this);
  double fb=(*eq_function)(zone_index,b,this);
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
    fb=(*eq_function)(zone_index,b,this);
  }
  printf("Maximum number of iterations exceeded in zbrent\n");
  return 0.0;
}
