/*
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

using namespace std;
namespace pc = physical_constants;

LorentzHelper::LorentzHelper(const double v_in[3], const bool exp_dec){
	exponential_decay = exp_dec;
	// poison the frame-dependent quantities
	for(int i=1; i<2; i++){
		dist[i] = -1;
		absopac[i] = -1;
		scatopac[i] = -1;
		p[i] = Particle();
	}

	// set the relativistic quantities to be used later
	for(int i=0; i<3; i++) v[i] = v_in[i];
}

//==========//
// particle //
//==========//

template<Frame f>
void LorentzHelper::set_p(const Particle* p_in){
	p[lab] = *p_in;
	p[com] = *p_in;

	double vrel[3]; // v_other - v_f
	for(int i=0; i<3; i++) vrel[i] = v[i] * (f==lab ? 1.0 : -1.0);

	Frame other = ( f==com ? lab : com );
	lorentz_transform_particle(&(p[other]), vrel, 3);
}
template void LorentzHelper::set_p<com>(const Particle* p);
template void LorentzHelper::set_p<lab>(const Particle* p);

void LorentzHelper::scale_p_e(const double factor){
	PRINT_ASSERT(factor,>=,0);
	p[lab].e *= factor;
	p[com].e *= factor;
}

void LorentzHelper::set_p_tau(const double tau){
	PRINT_ASSERT(tau,>=,0);
	p[lab].tau = tau;
	p[com].tau = tau;
}

Particle LorentzHelper::particle_copy(const Frame f) const{ return p[f];}
const Particle* LorentzHelper::particle_readonly(const Frame f) const { return &(p[f]);}
double LorentzHelper::p_e(const Frame f) const {return p[f].e;}
double LorentzHelper::p_nu(const Frame f) const {return p[f].nu;}
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
const double* LorentzHelper::p_x(const int size) const {
	PRINT_ASSERT(size,==,3);
	return p[lab].x;
}
const double* LorentzHelper::p_D(Frame f, const int size) const{
	PRINT_ASSERT(size,==,3);
	return p[f].D;
}
void LorentzHelper::set_p_x(const double x[3], const int size) {
	PRINT_ASSERT(size,==,3);
	for(int i=0; i<size; i++){
		p[com].x[i] = x[i];
		p[lab].x[i] = x[i];
	}
}
void LorentzHelper::set_p_fate(const ParticleFate fate){
	p[com].fate = fate;
	p[lab].fate = fate;
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
	double dshift = p[f].nu / p[other].nu;
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
	double dshift = p[f].nu / p[other].nu;
	dist[other] = dist[f] / dshift;
}
template void LorentzHelper::set_distance<com>(const double d);
template void LorentzHelper::set_distance<lab>(const double d);

//=======//
// other //
//=======//

double LorentzHelper::lorentz_factor(const double v[3], const int vsize) const{
	PRINT_ASSERT(vsize,<=,3);
	PRINT_ASSERT(Transport::dot(v,v,vsize),<,pc::c*pc::c);
	double beta2 = Transport::dot(v,v,vsize) * pc::inv_c * pc::inv_c;
	double lfac = 1.0 / sqrt(1.0 - beta2);
	PRINT_ASSERT(lfac,>=,1.0);
	return lfac;
}

double LorentzHelper::doppler_shift(const double v[3], const double D[3], const int size) const{
	// new frame is comoving frame. old frame is lab frame.
	// v_rel = v_comoving - v_lab  -->  v keeps its sign

	double gamma = lorentz_factor(v,size);
	double vdd = Transport::dot(v,D,size);
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
	PRINT_ASSERT(p->nu,>=,0);
	PRINT_ASSERT(p->e,>=,0);

	// calculate the doppler shift, v dot D, and lorentz factors
	double gamma = Transport::lorentz_factor(v,vsize);
	double vdd = Transport::dot(p->D, v, vsize);
	double dshift = doppler_shift(v, p->D, vsize);

	// transform the 0th component (energy and frequency)
	p->e  *= dshift;
	p->nu *= dshift;

	// transform the 1-3 components (direction)
	// See Mihalas & Mihalas eq 89.8
	double tmp = gamma/pc::c * (1 - gamma*vdd/(pc::c*(gamma+1)));
	p->D[0] = (p->D[0] - v[0]*tmp);
	p->D[1] = (p->D[1] - v[1]*tmp);
	p->D[2] = (p->D[2] - v[2]*tmp);
	Transport::normalize(p->D,3);

	// sanity checks
	PRINT_ASSERT(p->e,>=,0);
	PRINT_ASSERT(p->nu,>=,0);
	PRINT_ASSERT(dshift,>,0);
}
