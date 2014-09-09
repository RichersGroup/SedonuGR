#include <cassert>
#include "transport.h"
#include "particle.h"
#include "grid_general.h"
#include "global_options.h"

// lorentz factor ("gamma")
// v_rel = v_newframe - v_oldframe
double transport::lorentz_factor(const vector<double>& v){
	assert(v.size()==3);
	assert(dot(v,v) < pc::c*pc::c);
	double beta2 = dot(v,v) / (pc::c*pc::c);
	double lfac = 1.0 / sqrt(1.0 - beta2);
	assert(lfac>=1.0);
	return lfac;
}

// dot product of v_rel and relativistic particle direction
// v_rel = v_newframe - v_oldframe
// D = direction vector of relativistic particle in old frame
double transport::dot(const vector<double>& a, const vector<double>& b){
	assert(a.size()>0);
	assert(b.size()>0);
	assert(a.size()==b.size());
	double product = 0;
	for(unsigned i=0; i<a.size(); i++) product += a[i]*b[i];
	return product;
}

// v_dot_d is the dot product of the relative velocity and the relativistic particle's direction
double doppler_shift(const double gamma, const double vdd){
	assert(gamma > 0);
	double dshift = gamma * (1.0 - vdd/pc::c);
	assert(dshift>0);
	return dshift;
}

// apply a lorentz transform to the particle
// v = v_newframe - v_oldframe
void lorentz_transform(particle* p, const vector<double> v){
	// check input
	assert(v.size()==3);
	assert(p->nu > 0);
	assert(p->e > 0);

	// calculate the doppler shift, v dot D, and lorentz factors
	double gamma = transport::lorentz_factor(v);
	double vdd = transport::dot(v, p->D);
	double dshift = doppler_shift(gamma, vdd);

	// transform the 0th component (energy and frequency)
	p->e  *= dshift;
	p->nu *= dshift;

	// transform the 1-3 components (direction)
	// See Mihalas & Mihalas eq 89.8
	p->D[0] = 1.0/dshift * (p->D[0] - gamma*v[0]/pc::c * (1 - gamma*vdd/pc::c/(gamma+1)) );
	p->D[1] = 1.0/dshift * (p->D[1] - gamma*v[1]/pc::c * (1 - gamma*vdd/pc::c/(gamma+1)) );
	p->D[2] = 1.0/dshift * (p->D[2] - gamma*v[2]/pc::c * (1 - gamma*vdd/pc::c/(gamma+1)) );
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
	double vdd = dot(v, p->D);
	double dshift = doppler_shift(gamma,vdd);
	assert(dshift>0);
	return dshift;
}

double transport::dshift_lab_to_comoving(const particle* p) const
{
	vector<double> v;
	grid->cartesian_velocity_vector(p->x,v); // v_comoving - v_lab

	// new frame is comoving frame. old frame is lab frame.
	// v_rel = v_comoving - v_lab  -->  v keeps its sign

	double gamma = lorentz_factor(v);
	double vdd = dot(v, p->D);
	double dshift = doppler_shift(gamma,vdd);
	assert(dshift>0);
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

double transport::comoving_dt(const double lab_dt, const int z_ind) const{
	assert(lab_dt>0);
	assert(z_ind >= 0);
	assert(z_ind < (int)grid->z.size());
	return lab_dt / lorentz_factor(grid->z[z_ind].v);
}
