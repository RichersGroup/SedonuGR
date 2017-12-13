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

#include "global_options.h"
#include "LorentzHelper.h"
#include "Transport.h"
#include "physical_constants.h"
#include "Grid.h"

using namespace std;
namespace pc = physical_constants;

// initialize with the exponential_decay parameter
// so it knows which opacity to use for tau calculations
LorentzHelper::LorentzHelper(const bool exp_dec){
	exponential_decay = exp_dec;
	for(int i=1; i<2; i++){
		dist[i] = -1;
		absopac[i] = -1;
		scatopac[i] = -1;
		p[i] = Particle();
	}
	for(int i=0; i<3; i++) v[i] = 0;
}

// updates the fluid velocity, assuming the particle in the lab frame remains constant
void LorentzHelper::set_v(const double v_in[3], const int size){
	PRINT_ASSERT(size,==,3);

	// set the relativistic quantities to be used later
	for(int i=0; i<3; i++) v[i] = v_in[i];

	if(p[lab].N>=0 && p[lab].kup[3]>=0) set_p<lab>(&(p[lab]));
}

// get the velocity vector back out
const double* LorentzHelper::velocity(const int size) const{
	return v;
}

//==========//
// particle //
//==========//

// two ways to get a copy of the paricle. The readonly is more efficient.
Particle LorentzHelper::particle_copy(const Frame f) const{ return p[f];}
const Particle* LorentzHelper::particle_readonly(const Frame f) const { return &(p[f]);}


// set the particle quantities based on an input particle.
// Resets the distance and opacities if they have already been set
template<Frame f>
void LorentzHelper::set_p(const Particle* p_in){
	p[lab] = *p_in;
	p[com] = *p_in;

	double vrel[3]; // v_other - v_f
	for(int i=0; i<3; i++) vrel[i] = v[i] * (f==lab ? 1.0 : -1.0);

	Frame other = ( f==com ? lab : com );
	lorentz_transform_particle(&(p[other]), vrel, 3);

	// set the other things
	if(absopac[f]>=0 && scatopac[f]>=0) set_opac<f>(absopac[f],scatopac[f]);
	if(dist[f]>=0) set_distance<f>(dist[f]);
}
template void LorentzHelper::set_p<com>(const Particle* p);
template void LorentzHelper::set_p<lab>(const Particle* p);


// set the particle frame-dependent quantities individually
template<Frame f>
void LorentzHelper::set_p_kup(const double kup[4], const int size) {
	PRINT_ASSERT(size,==,4);
	for(int i=0; i<size; i++) p[f].kup[i] = kup[i];
	if(p[f].N>0) set_p<f>(&(p[f]));
}
template void LorentzHelper::set_p_kup<com>(const double kup[4], const int size);
template void LorentzHelper::set_p_kup<lab>(const double kup[4], const int size);


// rescale the particle energies. Does not require performing transformations
void LorentzHelper::scale_p_number(const double factor){
	PRINT_ASSERT(factor,>=,0);
	p[lab].N *= factor;
	p[com].N *= factor;
}
void LorentzHelper::scale_p_energy(const double factor){
	PRINT_ASSERT(factor,>=,0);
	for(int i=0; i<4; i++){
		p[lab].kup[i] *= factor;
		p[com].kup[i] *= factor;
	}
	scale_p_number(factor);
}


// set the lorentz-invariant particle properties
void LorentzHelper::set_p_tau(const double tau){
	PRINT_ASSERT(tau,>=,0);
	p[lab].tau = tau;
	p[com].tau = tau;
}
void LorentzHelper::set_p_xup(const double x[4], const int size) {
	PRINT_ASSERT(size,==,4);
	for(int i=0; i<size; i++){
		p[com].xup[i] = x[i];
		p[lab].xup[i] = x[i];
	}
}
void LorentzHelper::set_p_fate(const ParticleFate fate){
	p[com].fate = fate;
	p[lab].fate = fate;
}
void LorentzHelper::set_p_s(const int s){
	PRINT_ASSERT(s,>=,0);
	p[com].s = s;
	p[lab].s = s;
}


// get the particle properties
double LorentzHelper::p_N() const {
	PRINT_ASSERT(p[com].N,==,p[lab].N);
	return p[com].N;
}
double LorentzHelper::p_nu(const Frame f) const {return p[f].kup[3] * pc::c / (2.0*pc::pi);}
int    LorentzHelper::p_s() const {
	PRINT_ASSERT(p[lab].s,==,p[com].s);
	return p[lab].s;
}
double LorentzHelper::p_tau() const {
	PRINT_ASSERT(p[lab].tau,==,p[com].tau);
	return p[lab].tau;
}
ParticleFate LorentzHelper::p_fate() const{
	PRINT_ASSERT(p[lab].fate,==,p[com].fate);
	return p[lab].fate;
}
const double* LorentzHelper::p_xup() const {
	return p[lab].xup;
}
void LorentzHelper::p_D(const Frame f, double D[3], const int size) const{
	PRINT_ASSERT(size,==,3);
	PRINT_ASSERT(p[f].kup[3],>,0);
	for(int i=0; i<3; i++) D[i] = p[f].kup[i];
	Grid::normalize_Minkowski<3>(D,3);
}
const double* LorentzHelper::p_kup(const Frame f) const {
	return p[f].kup;
}



//=========//
// opacity //
//=========//

template<Frame f>
void LorentzHelper::set_opac(const double abs, const double scat){
	PRINT_ASSERT(abs,>=,0);
	PRINT_ASSERT(scat,>=,0);
	absopac[f] = abs;
	scatopac[f] = scat;

	// nu * opac is lorentz invariant
	Frame other = ( f==com ? lab : com );
	double dshift = p[f].kup[3] / p[other].kup[3];
	absopac[other]  =  absopac[f] * dshift;
	scatopac[other] = scatopac[f] * dshift;
}
template void LorentzHelper::set_opac<com>(const double abs, const double scat);
template void LorentzHelper::set_opac<lab>(const double abs, const double scat);

double LorentzHelper::net_opac (const Frame f) const {return scatopac[f] + absopac[f];}
double LorentzHelper::abs_opac (const Frame f) const {return  absopac[f]             ;}
double LorentzHelper::scat_opac(const Frame f) const {return scatopac[f]             ;}
double LorentzHelper::tau_opac(const Frame f) const {
	return exponential_decay ? scat_opac(f) : net_opac(f);
}
double LorentzHelper::abs_fraction() const {
	if(net_opac(com)==0) return 0;
	else return abs_opac(com) / net_opac(com);
}


//==========//
// distance //
//==========//

double LorentzHelper::distance(const Frame f) const {return dist[f];}

template<Frame f>
void LorentzHelper::set_distance(const double d){
	PRINT_ASSERT(d,>=,0);
	dist[f] = d;

	// l/nu is lorentz invariant (transforms oppositely as opac)
	Frame other = ( f==com ? lab : com );
	double dshift = p[f].kup[3] / p[other].kup[3];
	dist[other] = dist[f] / dshift;
}
template void LorentzHelper::set_distance<com>(const double d);
template void LorentzHelper::set_distance<lab>(const double d);

//=======//
// other //
//=======//




// apply a general lorentz transform to a 3D vector.
// first three components are spatial, 4th component is time (units of distance)
// input velocity is the fluid velocity in the lab frame
// [v] = cm/s, [x] = cm
void LorentzHelper::transform_cartesian_4vector_c2l(const double v[3], double x[4]){
	PRINT_ASSERT(Grid::dot_Minkowski<3>(v,v,3),<=,pc::c*pc::c);

	// new frame is lab frame. old frame is comoving frame.
	// v_rel = v_lab - v_comoving  --> v must flip sign.
	double vrel[3];
	vrel[0] = -v[0];
	vrel[1] = -v[1];
	vrel[2] = -v[2];

	double gamma = LorentzHelper::lorentz_factor(vrel,3);
	double v2 = Grid::dot_Minkowski<3>(vrel,vrel,3);


	// save comoving x
	double xcom[4];
	xcom[0] = x[0];
	xcom[1] = x[1];
	xcom[2] = x[2];
	xcom[3] = x[3];

	// time component of lab x
	x[4] = gamma*xcom[3];
	for(int i=0; i<3; i++) x[3] -= gamma*vrel[i]/pc::c * xcom[i];

	// spatial components
	for(int i=0; i<3; i++){
		x[i] = xcom[i];
		if(v2 > 0){
			x[i] -= gamma*vrel[i]/pc::c * xcom[3];
			for(int j=0; j<3; j++){
				x[i] += (gamma-1.0)*vrel[i]*vrel[j]/v2 * xcom[j];
			}
		}
	}

}

double LorentzHelper::lorentz_factor(const double v[3], const int vsize){
	PRINT_ASSERT(vsize,<=,3);
	PRINT_ASSERT(Grid::dot_Minkowski<3>(v,v,vsize),<,pc::c*pc::c);
	double beta2 = Grid::dot_Minkowski<3>(v,v,vsize) * pc::inv_c * pc::inv_c;
	double lfac = 1.0 / sqrt(1.0 - beta2);
	PRINT_ASSERT(lfac,>=,1.0);
	return lfac;
}

double LorentzHelper::doppler_shift(const double v[3], const double D[3], const int size) const{
	// new frame is comoving frame. old frame is lab frame.
	// v_rel = v_comoving - v_lab  -->  v keeps its sign

	double gamma = lorentz_factor(v,size);
	double vdd = Grid::dot_Minkowski<3>(v,D,size);
	double dshift = gamma * (1.0 - vdd*pc::inv_c);
	//double dshift = doppler_shift(gamma,vdd);
	PRINT_ASSERT(dshift,>,0);
	return dshift;
}


// apply a lorentz transform to the particle
// v = v_newframe - v_oldframe
// optimized for a null particle
void LorentzHelper::lorentz_transform_particle(Particle* p, const double v[3], const int vsize) const{
	// check input
	PRINT_ASSERT(vsize,==,3);
	PRINT_ASSERT(p->kup[3],>,0);
	PRINT_ASSERT(p->kup[3],<,INFINITY);
	PRINT_ASSERT(p->N,>=,0);

	// calculate the doppler shift, v dot D, and lorentz factors
	double D[3];
	for(int i=0; i<3; i++){
		PRINT_ASSERT(p->kup[i],<,INFINITY);
		D[i] = p->kup[i];
	}
	Grid::normalize_Minkowski<3>(D,3);
	double gamma = lorentz_factor(v,vsize);
	double vdd = Grid::dot_Minkowski<3>(D, v, vsize);
	PRINT_ASSERT(D[0],<,INFINITY);
	double dshift = doppler_shift(v, D, vsize);

	// transform the 0th component (frequency)
	p->kup[3] *= dshift;

	// transform the 1-3 components (direction)
	// See Mihalas & Mihalas eq 89.8
	double tmp = gamma/pc::c * (1 - gamma*vdd/(pc::c*(gamma+1)));
	D[0] = (D[0] - v[0]*tmp);
	D[1] = (D[1] - v[1]*tmp);
	D[2] = (D[2] - v[2]*tmp);
	Grid::normalize_Minkowski<3>(D,3);
	for(int i=0; i<3; i++) p->kup[i] = D[i] * p->kup[3];

	// sanity checks
	PRINT_ASSERT(p->N,>=,0);
	PRINT_ASSERT(p->kup[3],>=,0);
	PRINT_ASSERT(dshift,>,0);
}
