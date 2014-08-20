#include <cassert>
#include "transport.h"
#include "particle.h"
#include "physical_constants.h"
#include "grid_general.h"

namespace pc = physical_constants;

// lorentz factor ("gamma")
// v_rel = v_newframe - v_oldframe
double lorentz_factor(const vector<double> v_rel){
	assert(v_rel.size()==3);
  double beta2 = (v_rel[0]*v_rel[0] + v_rel[1]*v_rel[1] + v_rel[2]*v_rel[2]) / (pc::c*pc::c);
  double lfac = 1.0 / sqrt(1.0 - beta2);
  assert(lfac>=1.0);
  return lfac;
}

// dot product of v_rel and relativistic particle direction
// v_rel = v_newframe - v_oldframe
// D = direction vector of relativistic particle in old frame
double v_dot_d(const vector<double> v_rel, const vector<double> D){
	assert(v_rel.size()==3);
	assert(D.size()==3);
  return v_rel[0]*D[0] + v_rel[1]*D[1] + v_rel[2]*D[2];
}

// v_dot_d is the dot product of the relative velocity and the relativistic particle's direction
double doppler_shift(const double gamma, const double vdd){
	assert(gamma > 0);
  double dshift = gamma * (1.0 - vdd/pc::c);
  assert(dshift>0);
  return dshift;
}

// apply a lorentz transform to the particle
// v_rel = v_newframe - v_oldframe
void lorentz_transform(particle* p, const vector<double> v_rel){
  // check input
	assert(v_rel.size()==3);
  assert(p->nu > 0);
  assert(p->e > 0);

  // calculate the doppler shift, v dot D, and lorentz factors
  double gamma = lorentz_factor(v_rel);
  double vdd = v_dot_d(v_rel, p->D);
  double dshift = doppler_shift(gamma, vdd);

  // transform the 0th component (energy and frequency)
  p->e  *= dshift;
  p->nu *= dshift;

  // transform the 1-3 components (direction)
  // See Mihalas & Mihalas eq 89.8
  p->D[0] = 1.0/dshift * (p->D[0] - gamma*v_rel[0]/pc::c * (1 - gamma*vdd/pc::c/(gamma+1)) );
  p->D[1] = 1.0/dshift * (p->D[1] - gamma*v_rel[1]/pc::c * (1 - gamma*vdd/pc::c/(gamma+1)) );
  p->D[2] = 1.0/dshift * (p->D[2] - gamma*v_rel[2]/pc::c * (1 - gamma*vdd/pc::c/(gamma+1)) );
  p->normalize_direction();

  // sanity checks
  assert(p->e > 0);
  assert(p->nu > 0);
  assert(dshift > 0);
}


//------------------------------------------------------------
// get the doppler shift when moving from frame_to_frame
// does not change any particle properties
//------------------------------------------------------------
double transport::dshift_comoving_to_lab(const particle* p) const
{
  vector<double> v;
  grid->cartesian_velocity_vector(p->x,v); // v_comoving - v_lab

  // new frame is lab frame. old frame is comoving frame.
  // v_rel = v_lab - v_comoving  --> v must flip sign.
  v[0] *= -1;
  v[1] *= -1;
  v[2] *= -1;

  double gamma = lorentz_factor(v);
  double vdd = v_dot_d(v, p->D);
  double dshift = doppler_shift(gamma,vdd);
  return dshift;
}

double transport::dshift_lab_to_comoving(const particle* p) const
{
	vector<double> v;
  grid->cartesian_velocity_vector(p->x,v); // v_comoving - v_lab

  // new frame is comoving frame. old frame is lab frame.
  // v_rel = v_comoving - v_lab  -->  v keeps its sign

  double gamma = lorentz_factor(v);
  double vdd = v_dot_d(v, p->D);
  double dshift = doppler_shift(gamma,vdd);
  return dshift;
}


//------------------------------------------------------------
// do a lorentz transformation; modifies the energy, frequency
// and direction vector of the particle
//------------------------------------------------------------
void transport::transform_comoving_to_lab(particle* p) const
{
  vector<double> v;
  grid->cartesian_velocity_vector(p->x,v); // v_comoving - v_lab

  // new frame is lab frame. old frame is comoving frame.
  // v_rel = v_lab - v_comoving  --> v must flip sign.
  v[0] *= -1;
  v[1] *= -1;
  v[2] *= -1;

  lorentz_transform(p,v);
}

void transport::transform_lab_to_comoving(particle* p) const
{
  vector<double> v;
  grid->cartesian_velocity_vector(p->x,v); // v_comoving - v_lab

  // new frame is lab frame. old frame is comoving frame.
  // v_rel = v_comoving - v_lab  --> v keeps its sign.

  lorentz_transform(p,v);
}
