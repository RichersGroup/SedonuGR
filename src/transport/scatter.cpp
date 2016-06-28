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
void Transport::event_interact(Particle* p, const int z_ind, const double abs_frac, const double lab_opac, const double com_opac){
	PRINT_ASSERT(z_ind,>=,0);
	PRINT_ASSERT(z_ind,<,(int)grid->z.size());
	PRINT_ASSERT(abs_frac,>=,0.0);
	PRINT_ASSERT(abs_frac,<=,1.0);
	PRINT_ASSERT(p->e,>,0);

	int do_absorb_partial  = !radiative_eq;
	int do_absorb_reemit   =  radiative_eq;

	// particle is transformed to the comoving frame
	transform_lab_to_comoving(p,z_ind);

	// absorb part of the packet
	if(do_absorb_partial){
		p->e *= (1.0 - abs_frac);
		window(p,z_ind);
		if(p->fate==moving) scatter(p,abs_frac,com_opac, z_ind);
	}
	// absorb the particle and let the fluid re-emit another particle
	else if(do_absorb_reemit){
		re_emit(p,z_ind);
		L_net_lab[p->s] += p->e;
		PRINT_ASSERT(p->e,>,0.0);
	}

	if(p->fate==moving){
		// particle is transformed back to the lab frame
		PRINT_ASSERT(p->e,>,0);
		transform_comoving_to_lab(p,z_ind);

		// resample the path length
		sample_tau(p,lab_opac,abs_frac);
		if(p->fate==moving) window(p,z_ind);

		// sanity checks
		if(p->fate==moving){
			PRINT_ASSERT(p->nu,>,0);
			PRINT_ASSERT(p->e,>,0);
		}
	}
}


// decide whether to kill a particle
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
void Transport::re_emit(Particle* p, const int z_ind) const{
	PRINT_ASSERT(z_ind,>=,0);
	PRINT_ASSERT(z_ind,<,(int)grid->z.size());

	// reset the particle properties
	isotropic_direction(p);
	sample_zone_species(p,z_ind);
	species_list[p->s]->sample_zone_nu(*p,z_ind);

	// tally into zone's emitted energy
	grid->z[z_ind].e_emit += p->e;

	// sanity checks
	PRINT_ASSERT(p->nu,>,0);
	PRINT_ASSERT(p->s,>=,0);
	PRINT_ASSERT(p->s,<,(int)species_list.size());
}

// choose which type of scattering event to do
void Transport::scatter(Particle* p_comoving, double abs_frac, double com_opac, int z_ind) const{
	PRINT_ASSERT(abs_frac,>=,0);
	PRINT_ASSERT(abs_frac,<=,1.0);
	PRINT_ASSERT(com_opac,>=,0);
	bool did_random_walk = false;

	// try to do random walk approximation in scattering-dominated diffusion regime
	if(randomwalk_sphere_size>0){

		double com_absopac = com_opac * abs_frac;
		double com_scatopac = com_opac * (1.0-abs_frac);
		double D = pc::c / (3.0 * com_scatopac); // diffusion coefficient (cm^2/s)

		// if the optical depth is below our threshold, don't do random walk
		// (first pass to avoid doing lots of math)
		double Rlab = randomwalk_sphere_size * grid->zone_min_length(z_ind);
		if(com_scatopac * Rlab >= randomwalk_min_optical_depth){
			// determine maximum comoving sphere size
			double v[3] = {0,0,0};
			grid->cartesian_velocity_vector(p_comoving->x,3,v,3,z_ind);
			double vabs = sqrt(dot(v,v,3));
			double gamma = lorentz_factor(v,3);

			double Rcom = 0;
			if(Rlab==0) Rcom = 0;
			else if(Rlab==INFINITY) Rcom = randomwalk_min_optical_depth / (com_absopac>0 ? com_absopac : com_scatopac);
			else Rcom =  2. * Rlab / gamma / (1. + sqrt(1. + 4.*Rlab*vabs*randomwalk_max_x / (gamma*D) ) );

			// if the optical depth is below our threshold, don't do random walk
			if(com_scatopac * Rcom >= randomwalk_min_optical_depth){
				double eold = p_comoving->e;
				random_walk(p_comoving,com_absopac,com_scatopac, Rcom, D, z_ind);
				did_random_walk = true;
			}
		}
	}

	// isotropic scatter if can't do random walk
	if(!did_random_walk) isotropic_direction(p_comoving);

}


// isotropic scatter, done in COMOVING frame
void Transport::isotropic_direction(Particle* p_comoving) const
{
	// Randomly generate new direction isotropically in comoving frame
	double mu  = 1 - 2.0*rangen.uniform();
	double phi = 2.0*pc::pi*rangen.uniform();
	double smu = sqrt(1 - mu*mu);
	p_comoving->D[0] = smu*cos(phi);
	p_comoving->D[1] = smu*sin(phi);
	p_comoving->D[2] = mu;
	normalize(p_comoving->D,3);
}

//---------------------------------------------------------------------
// Randomly select an optical depth through which a particle will move.
//---------------------------------------------------------------------
void Transport::sample_tau(Particle *p, const double lab_opac, const double abs_frac){
	PRINT_ASSERT(bias_path_length,>=,0);
	PRINT_ASSERT(p->fate,==,moving);

	// tweak distribution if biasing path length
	// this is probably computed more times than necessary. OPTIMIZE
	double taubar = 1.0;
	if(bias_path_length && abs_frac<1.0){
		taubar = 1.0/(1.0-abs_frac);
		taubar = min(taubar,max_path_length_boost);
	}

	// sample the distribution and modify the energy
	do{ // don't allow tau to be infinity
		p->tau = -taubar*log(rangen.uniform());
	} while(p->tau >= INFINITY);
	if(taubar != 1.0) p->e *= taubar * exp(-p->tau * (1.0 - 1.0/taubar));
	if(p->e==0) p->fate = rouletted;
	PRINT_ASSERT(p->e,>=,0);

	// make sure nothing crazy happened
	if(p->fate==moving){
		PRINT_ASSERT(p->tau,>=,0);
		PRINT_ASSERT(p->e,>,0);
		PRINT_ASSERT(p->e,<,INFINITY);
		PRINT_ASSERT(p->tau,<,INFINITY);
	}
	else PRINT_ASSERT(p->fate,==,rouletted);
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
	double displacement4[4] = {Rcom*mu*cos(phi), Rcom*mu*sin(phi), Rcom*(1.0-mu*mu), distance};


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
	transform_cartesian_4vector_c2l(grid->z[z_ind].u, displacement4);
	double displacement3[3] = {displacement4[0], displacement4[1], displacement4[2]};
	normalize(displacement3,3);
	double dshift = dshift_comoving_to_lab(p->x,displacement3,z_ind);
	zone->distribution[p->s].count(displacement3, 3, p->nu*dshift, erad_com*dshift);

	// move the particle to the edge of the sphere and transform to the lab frame
	for(int i=0; i<3; i++) p->x[i] += displacement3[3];
	if(ratio_deposited > 0) p->e *= 1.0 - ratio_deposited;
	if(p->e == 0) p->fate = rouletted;
}
