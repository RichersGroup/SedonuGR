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

#include "Species.h"
#include "Transport.h"
#include "Grid.h"
#include "global_options.h"

using namespace std;
namespace pc = physical_constants;

//------------------------------------------------------------
// physics of absorption/scattering
//------------------------------------------------------------
void Transport::event_interact(LorentzHelper* lh, const int z_ind){
	PRINT_ASSERT(z_ind,>=,0);
	PRINT_ASSERT(z_ind,<,(int)grid->z.size());
	PRINT_ASSERT(lh->abs_fraction(),>=,0.0);
	PRINT_ASSERT(lh->abs_fraction(),<=,1.0);
	PRINT_ASSERT(lh->p_e(com),>,0);
	PRINT_ASSERT(lh->p_fate(),==,moving);

	// absorb the particle and let the fluid re-emit another particle
	if(radiative_eq){
		re_emit(lh,z_ind);
		L_net_lab[lh->p_s()] += lh->p_e(lab);
	}

	// absorb part of the packet
	else{
		if(!exponential_decay) lh->scale_p_e(1.0 - lh->abs_fraction());
		scatter(lh, z_ind);
	}

	// resample the path length
	if(lh->p_fate()==moving) sample_tau(lh);

	// window the particle
	if(lh->p_fate()==moving) window(lh,z_ind);

	// sanity checks
	if(lh->p_fate()==moving){
		PRINT_ASSERT(lh->p_nu(com),>,0);
		PRINT_ASSERT(lh->p_nu(com),>,0);
	}
}


// decide whether to kill a particle
void Transport::window(LorentzHelper *lh, const int z_ind){
	PRINT_ASSERT(lh->p_e(com),>=,0);
	PRINT_ASSERT(lh->p_fate(),!=,rouletted);

	// Roulette if too low energy
	while(lh->p_e(com)<=min_packet_energy && lh->p_fate()==moving){
		if(rangen.uniform() < 0.5) lh->set_p_fate(rouletted);
		else lh->scale_p_e(2.0);
	}
	if(lh->p_fate()==moving) PRINT_ASSERT(lh->p_e(com),>=,min_packet_energy);

	// split if too high energy, if enough space, and if in important region
	double ratio = lh->p_e(com) / max_packet_energy;
	int n_new = (int)ratio;
	if(ratio>1.0 && particles.size()+n_new<max_particles && species_list[lh->p_s()]->interpolate_importance(lh->p_nu(com),z_ind)>=1.0){
		lh->scale_p_e( 1.0 / (double)(n_new+1) );
		for(int i=0; i<n_new; i++){
			#pragma omp critical
			particles.push_back(lh->particle_copy(lab));
		}
	}

	if(lh->p_fate() == moving){
		PRINT_ASSERT(lh->p_e(com),<,INFINITY);
		PRINT_ASSERT(lh->p_e(com),>,0);
	}
	if(particles.size()>=max_particles && verbose && rank0){
		cout << "max_particles: " << max_particles << endl;
		cout << "particles.size(): " << particles.size() << endl;
		cout << "WARNING: max_particles is too small to allow splitting." << endl;
	}
}
void Transport::window(Particle* p, const int z_ind){
	PRINT_ASSERT(p->fate,!=,rouletted);

	// Roulette if too low energy
	while(p->e<=min_packet_energy && p->fate==moving){
		if(rangen.uniform() < 0.5) p->fate = rouletted;
		else p->e *= 2.0;
	}
	if(p->fate==moving) PRINT_ASSERT(p->e,>=,min_packet_energy);

	// split if too high energy, if enough space, and if in important region
	while(p->e>max_packet_energy && particles.size()<max_particles && species_list[p->s]->interpolate_importance(p->nu,z_ind)>=1.0){
		p->e /= 2.0;
		Particle pnew = *p;
		window(&pnew,z_ind);
		#pragma omp critical
		particles.push_back(pnew);
	}
	if(p->fate == moving){
		PRINT_ASSERT(p->e,<,INFINITY);
		PRINT_ASSERT(p->e,>,0);
	}
	if(particles.size()>=max_particles && verbose && rank0){
		cout << "max_particles: " << max_particles << endl;
		cout << "particles.size(): " << particles.size() << endl;
		cout << "WARNING: max_particles is too small to allow splitting." << endl;
	}
}

// re-emission, done in COMOVING frame
void Transport::re_emit(LorentzHelper *lh, const int z_ind) const{
	PRINT_ASSERT(z_ind,>=,0);
	PRINT_ASSERT(z_ind,<,(int)grid->z.size());

	// reset the particle properties
	double Dnew[3];
	isotropic_direction(Dnew,3);
	lh->set_p_D<com>(Dnew,3);
	sample_zone_species(lh,z_ind);
	species_list[lh->p_s()]->sample_zone_nu(lh,z_ind);

	// tally into zone's emitted energy
	grid->z[z_ind].e_emit += lh->p_e(com);

	// sanity checks
	PRINT_ASSERT(lh->p_nu(com),>,0);
	PRINT_ASSERT(lh->p_s(),>=,0);
	PRINT_ASSERT(lh->p_s(),<,(int)species_list.size());
}

// choose which type of scattering event to do
void Transport::scatter(LorentzHelper *lh, int z_ind) const{\
	PRINT_ASSERT(lh->abs_fraction(),>=,0);
	PRINT_ASSERT(lh->abs_fraction(),<=,1.0);
	PRINT_ASSERT(lh->net_opac(com),>=,0);
	bool did_random_walk = false;

	// try to do random walk approximation in scattering-dominated diffusion regime
	if(randomwalk_sphere_size>0){

		double D = pc::c / (3.0 * lh->scat_opac(com)); // diffusion coefficient (cm^2/s)

		// if the optical depth is below our threshold, don't do random walk
		// (first pass to avoid doing lots of math)
		double Rlab = randomwalk_sphere_size * grid->zone_min_length(z_ind);
		if(lh->scat_opac(com) * Rlab >= randomwalk_min_optical_depth){
			// determine maximum comoving sphere size
			const double* v = lh->velocity(3);
			double vabs = sqrt(dot(v,v,3));
			double gamma = lorentz_factor(v,3);

			double Rcom = 0;
			if(Rlab==0) Rcom = 0;
			else if(Rlab==INFINITY) Rcom = randomwalk_min_optical_depth / (lh->abs_opac(com)>0 ? lh->abs_opac(com) : lh->scat_opac(com));
			else Rcom =  2. * Rlab / gamma / (1. + sqrt(1. + 4.*Rlab*vabs*randomwalk_max_x / (gamma*D) ) );

			// if the optical depth is below our threshold, don't do random walk
			if(lh->scat_opac(com) * Rcom >= randomwalk_min_optical_depth){
				{
				Particle p = lh->particle_copy(com);
				random_walk(&p,lh->abs_opac(com),lh->scat_opac(com), Rcom, D, z_ind);
				lh->set_p<com>(&p);
				}
				did_random_walk = true;
			}
		}
	}

	// isotropic scatter if can't do random walk
	if(!did_random_walk){
		double Dnew[3];
		isotropic_direction(Dnew,3);
		lh->set_p_D<com>(Dnew,3);
	}
}


// isotropic scatter, done in COMOVING frame
void Transport::isotropic_direction(double D[3], const int size) const
{
	PRINT_ASSERT(size,==,3);

	// Randomly generate new direction isotropically in comoving frame
	double mu  = 1 - 2.0*rangen.uniform();
	double phi = 2.0*pc::pi*rangen.uniform();
	double smu = sqrt(1 - mu*mu);

	D[0] = smu*cos(phi);
	D[1] = smu*sin(phi);
	D[2] = mu;
	normalize(D,3);
}

//---------------------------------------------------------------------
// Randomly select an optical depth through which a particle will move.
//---------------------------------------------------------------------
void Transport::sample_tau(LorentzHelper *lh){
	PRINT_ASSERT(bias_path_length,>=,0);
	PRINT_ASSERT(lh->p_fate(),==,moving);

	// tweak distribution if biasing path length
	// this is probably computed more times than necessary. OPTIMIZE
	double taubar = 1.0;
	if(bias_path_length && lh->abs_fraction()<1.0){
		taubar = 1.0/(1.0-lh->abs_fraction());
		taubar = min(taubar,max_path_length_boost);
	}

	// sample the distribution and modify the energy
	do{ // don't allow tau to be infinity
		lh->set_p_tau( -taubar*log(rangen.uniform()) );
	} while(lh->p_tau() >= INFINITY);
	if(taubar != 1.0) lh->scale_p_e( taubar * exp(-lh->p_tau() * (1.0 - 1.0/taubar)) );
	if(lh->p_e(lab)==0) lh->set_p_fate(rouletted);
	PRINT_ASSERT(lh->p_e(lab),>=,0);

	// make sure nothing crazy happened
	if(lh->p_fate()==moving){
		PRINT_ASSERT(lh->p_tau(),>=,0);
		PRINT_ASSERT(lh->p_e(lab),>,0);
		PRINT_ASSERT(lh->p_e(lab),<,INFINITY);
		PRINT_ASSERT(lh->p_tau(),<,INFINITY);
	}
	else PRINT_ASSERT(lh->p_fate(),==,rouletted);
}


//-------------------------------------------------------
// Initialize the CDF that determines particle dwell time
// result is D*t/(R^2)
//-------------------------------------------------------
void Transport::init_randomwalk_cdf(Lua* lua){
	int sumN = lua->scalar<int>("randomwalk_sumN");
	int npoints = lua->scalar<int>("randomwalk_npoints");
	double max_x = lua->scalar<double>("randomwalk_max_x");
	double interpolation_order = lua->scalar<double>("randomwalk_interpolation_order");

	randomwalk_diffusion_time.resize(npoints);
	randomwalk_diffusion_time.interpolation_order = interpolation_order;
	randomwalk_xaxis.init(0,max_x,npoints, linear);

	#pragma omp parallel for
	for(int i=1; i<=npoints; i++){
		double sum = 0;
		double x = randomwalk_xaxis.x[i];

		for(int n=1; n<=sumN; n++){
			double tmp = 2.0/(double)n * exp(-x * (n*pc::pi)*(n*pc::pi));
			if(n%2 == 0) tmp *= -1;
			sum += tmp;
	    }

		randomwalk_diffusion_time.set(i,sum);
	}
	randomwalk_diffusion_time.normalize(0);
}

//----------------------
// Do a random walk step
//----------------------
void Transport::random_walk(Particle* p, const double com_absopac, const double com_scatopac, const double Rcom, const double D, const int z_ind) const{
	PRINT_ASSERT(com_scatopac,>,0);
	PRINT_ASSERT(com_absopac,>=,0);

	// set pointer to the current zone
	Zone* zone;
	zone = &(grid->z[z_ind]);


	// sample the distance travelled during the random walk
	double distance = pc::c * Rcom*Rcom / D * randomwalk_diffusion_time.invert(rangen.uniform(),&randomwalk_xaxis,-1);
	PRINT_ASSERT(distance,>=,Rcom);

	// deposit fluid quantities (comoving frame)
	double ratio_deposited = 1.0 - exp(-com_absopac * distance);
	#pragma omp atomic
	zone->e_abs += p->e * ratio_deposited;
	double this_l_comoving = 0;
	if(species_list[p->s]->lepton_number != 0){
		this_l_comoving = p->e/(p->nu*pc::h);
		double to_add = this_l_comoving * ratio_deposited;
		if(species_list[p->s]->lepton_number == 1){
            #pragma omp atomic
			zone->nue_abs += to_add;
		}
		else if(species_list[p->s]->lepton_number == -1){
            #pragma omp atomic
			zone->anue_abs += to_add;
		}
	}

	// randomly place the particle somewhere on the sphere (comoving frame)
	double phi = 2*pc::pi * rangen.uniform();
	double mu = 2.0*rangen.uniform() - 1.0;
	double displacement4[4] = {Rcom*mu*cos(phi), Rcom*mu*sin(phi), Rcom*sqrt(1.0-mu*mu), distance};


	//------------------------------------------------------------------------
	// pick a random outward direction, starting distribution pointing along z (comoving frame)
	phi = 2*pc::pi * rangen.uniform();
	mu = rangen.uniform();
	p->D[0] = mu*cos(phi);
	p->D[1] = mu*sin(phi);
	p->D[2] = sqrt(1.0-mu*mu);

	// get the displacement vector polar coordinates
	double d3com[3] = {displacement4[0], displacement4[1], displacement4[2]};
	normalize(d3com,3);
	double costheta = d3com[2];
	double sintheta = sqrt(d3com[0]*d3com[0] + d3com[1]*d3com[1]);

	if(abs(sintheta) < grid->tiny) p->D[2] *= costheta>0 ? 1.0 : -1.0;
	else{
		double cosphi = d3com[0] / sintheta;
		double sinphi = d3com[1] / sintheta;

		// first rotate away from the z axis along y=0 (move it toward x=0)
		normalize(p->D,3);
		double tmp[3];
		for(int i=0; i<3; i++) tmp[i] = p->D[i];
		p->D[0] =  costheta*tmp[0] + sintheta*tmp[2];
		p->D[2] = -sintheta*tmp[0] + costheta*tmp[2];

		// second rotate around the z axis, away from the x-axis
		normalize(p->D,3);
		for(int i=0; i<3; i++) tmp[i] = p->D[i];
		p->D[0] = cosphi*tmp[0] - sinphi*tmp[1];
		p->D[1] = sinphi*tmp[0] + cosphi*tmp[1];
	}
	normalize(p->D,3);
	//------------------------------------------------------------------------

	// calculate radiation energy in the comoving frame
	//double e_rad_directional = e_avg * R;
	//double e_rad_each_bin = e_avg * (distance - R) / (double)(zone->distribution[p->s].phi_dim() * zone->distribution[p->s].mu_dim());
	double e_avg = p->e / (com_absopac * distance) * (1.0 - exp(-com_absopac * distance));
	double erad_com = e_avg * distance;

	// deposit all radiaton energy into the bin corresponding to the direction of motion.
	// really, most should be isotropic and some should be in the direction of motion,
	// but this should average out properly over many trajectories.
	// depositing radiation in every bin would lead to lots of memory contention
	//normalize(d3lab,3);
	Particle pfake = *p;
	pfake.e = erad_com;
	for(int i=0; i<3; i++) pfake.D[i] = d3com[i];
	transform_comoving_to_lab(&pfake,z_ind);
	zone->distribution[p->s].count(pfake.D, 3, pfake.nu, pfake.e);

	// move the particle to the edge of the sphere
	transform_cartesian_4vector_c2l(zone->u, displacement4);
	double d3lab[3] = {displacement4[0], displacement4[1], displacement4[2]};
	for(int i=0; i<3; i++) p->x[i] += d3lab[3];
	if(ratio_deposited > 0) p->e *= 1.0 - ratio_deposited;
}
