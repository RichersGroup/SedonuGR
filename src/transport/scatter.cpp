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
	if(ratio>1.0 && (int)particles.size()+n_new<max_particles && species_list[lh->p_s()]->interpolate_importance(lh->p_nu(com),z_ind)>=1.0){
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
	if((int)particles.size()>=max_particles && verbose && rank0){
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
	double nu = p->kup[3]*pc::c/pc::h;
	while(p->e>max_packet_energy && (int)particles.size()<max_particles && species_list[p->s]->interpolate_importance(nu,z_ind)>=1.0){
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
	if((int)particles.size()>=max_particles && verbose && rank0){
		cout << "max_particles: " << max_particles << endl;
		cout << "particles.size(): " << particles.size() << endl;
		cout << "WARNING: max_particles is too small to allow splitting." << endl;
	}
}

// re-emission, done in COMOVING frame
void Transport::re_emit(LorentzHelper *lh, const int z_ind) const{
	PRINT_ASSERT(z_ind,>=,0);
	PRINT_ASSERT(z_ind,<,(int)grid->z.size());

	//double Ep = lh->p_e(com);
	Particle p = lh->particle_copy(com);

	// reset the particle properties
	p.s = sample_zone_species(z_ind,&p.e);
	double nu = species_list[p.s]->sample_zone_nu(z_ind,&p.e);
	grid->isotropic_kup(nu,p.kup,lh->p_xup(),4,&rangen);

	// set the LorentzHelper
	lh->set_p<com>(&p);

	// tally into zone's emitted energy
	grid->z[z_ind].e_emit += lh->p_e(com);

	// sanity checks
	PRINT_ASSERT(lh->p_nu(com),>,0);
	PRINT_ASSERT(lh->p_s(),>=,0);
	PRINT_ASSERT(lh->p_s(),<,(int)species_list.size());
}

// choose which type of scattering event to do
void Transport::scatter(LorentzHelper *lh, int z_ind) const{
	PRINT_ASSERT(lh->abs_fraction(),>=,0);
	PRINT_ASSERT(lh->abs_fraction(),<=,1.0);
	PRINT_ASSERT(lh->net_opac(com),>=,0);
	bool did_random_walk = false;

	// try to do random walk approximation in scattering-dominated diffusion regime
	if(randomwalk_sphere_size>0){

		double D = pc::c / (3.0 * lh->scat_opac(com)); // diffusion coefficient (cm^2/s)

		// if the optical depth is below our threshold, don't do random walk
		// (first pass to avoid doing lots of math)
		double Rlab_min = randomwalk_sphere_size * grid->zone_min_length(z_ind);
		double Rlab_boundary = grid->zone_cell_dist(lh->p_xup(),z_ind);
		double Rlab = max(Rlab_min, Rlab_boundary);
		if(lh->scat_opac(com) * Rlab >= randomwalk_min_optical_depth){
			// determine maximum comoving sphere size
			const double* v = lh->velocity(3);
			double vabs = sqrt(Grid::dot_Minkowski<3>(v,v,3));
			double gamma = LorentzHelper::lorentz_factor(v,3);

			double Rcom = 0;
			if(Rlab==0) Rcom = 0;
			else if(Rlab==INFINITY) Rcom = randomwalk_sphere_size * randomwalk_min_optical_depth / (lh->abs_opac(com)>0 ? lh->abs_opac(com) : lh->scat_opac(com));
			else Rcom =  2. * Rlab / gamma / (1. + sqrt(1. + 4.*Rlab*vabs*randomwalk_max_x / (gamma*D) ) );

			// if the optical depth is below our threshold, don't do random walk
			if(lh->scat_opac(com) * Rcom >= randomwalk_min_optical_depth){
				random_walk(lh, Rcom, D, z_ind);
				int new_ind = grid->zone_index(lh->p_xup(),3);
				if(new_ind==-2) lh->set_p_fate(escaped);
				did_random_walk = true;
			}
		}
	}

	// isotropic scatter if can't do random walk
	if(!did_random_walk && lh->p_fate()==moving){
		double kup[4];
		grid->isotropic_kup(lh->p_nu(com),kup,lh->p_xup(),4,&rangen);
		lh->set_p_kup<com>(kup,4);
	}
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
	randomwalk_max_x = lua->scalar<double>("randomwalk_max_x");
	double interpolation_order = lua->scalar<double>("randomwalk_interpolation_order");

	randomwalk_diffusion_time.resize(npoints);
	randomwalk_diffusion_time.interpolation_order = interpolation_order;
	randomwalk_xaxis.init(0,randomwalk_max_x,npoints, linear);

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
void Transport::random_walk(LorentzHelper *lh, const double Rcom, const double D, const int z_ind) const{
	PRINT_ASSERT(lh->scat_opac(com),>,0);
	PRINT_ASSERT(lh->abs_opac(com),>=,0);
	PRINT_ASSERT(lh->p_e(com),>=,0);
	PRINT_ASSERT(lh->p_nu(com),>=,0);

	// set pointer to the current zone
	Zone* zone;
	zone = &(grid->z[z_ind]);


	// sample the distance travelled during the random walk
	double path_length_com = pc::c * Rcom*Rcom / D * randomwalk_diffusion_time.invert(rangen.uniform(),&randomwalk_xaxis,-1);
	path_length_com = max(path_length_com,Rcom);

	// deposit fluid quantities (comoving frame)
	double ratio_deposited = 1.0 - exp(-lh->abs_opac(com) * path_length_com);
	#pragma omp atomic
	zone->e_abs += lh->p_e(com) * ratio_deposited;
	double this_l_comoving = 0;
	if(species_list[lh->p_s()]->lepton_number != 0){
		this_l_comoving = lh->p_e(com)/(lh->p_nu(com)*pc::h);
		double to_add = this_l_comoving * ratio_deposited;
		if(species_list[lh->p_s()]->lepton_number == 1){
            #pragma omp atomic
			zone->nue_abs += to_add;
		}
		else if(species_list[lh->p_s()]->lepton_number == -1){
            #pragma omp atomic
			zone->anue_abs += to_add;
		}
	}

	// randomly place the particle somewhere on the sphere (comoving frame)
	double Diso[3];
	grid->isotropic_direction(Diso,3,&rangen);
	double displacement4[4] = {Rcom*Diso[0], Rcom*Diso[1], Rcom*Diso[2], path_length_com};
	double d3com[3] = {displacement4[0], displacement4[1], displacement4[2]};
	LorentzHelper::transform_cartesian_4vector_c2l(zone->u, displacement4);
	double d3lab[3] = {displacement4[0], displacement4[1], displacement4[2]};

	//------------------------------------------------------------------------
	// pick a random outward direction, starting distribution pointing along z (comoving frame)
	double pD[3];
	double outward_phi = 2*pc::pi * rangen.uniform();
	double outward_mu = rangen.uniform();
	pD[0] = outward_mu*cos(outward_phi);
	pD[1] = outward_mu*sin(outward_phi);
	pD[2] = sqrt(1.0 - outward_mu*outward_mu);

	// get the displacement vector polar coordinates
	// theta_rotate is angle away from z-axis, not from xy-plane
	Grid::normalize_Minkowski<3>(d3com,3);
	double costheta_rotate = d3com[2];
	double sintheta_rotate = sqrt(1.0 - d3com[2]*d3com[2]);

	if(abs(sintheta_rotate) < TINY) pD[2] *= costheta_rotate>0 ? 1.0 : -1.0;
	else{

		// first rotate away from the z axis along y=0 (move it toward x=0)
		Grid::normalize_Minkowski<3>(pD,3);
		double pD_old[3];
		for(int i=0; i<3; i++) pD_old[i] = pD[i];
		pD[0] =  costheta_rotate*pD_old[0] + sintheta_rotate*pD_old[2];
		pD[2] = -sintheta_rotate*pD_old[0] + costheta_rotate*pD_old[2];

		// second rotate around the z axis, away from the x-axis
		double cosphi_rotate = d3com[0] / sintheta_rotate;
		double sinphi_rotate = d3com[1] / sintheta_rotate;
		Grid::normalize_Minkowski<3>(pD,3);
		for(int i=0; i<3; i++) pD_old[i] = pD[i];
		pD[0] = cosphi_rotate*pD_old[0] - sinphi_rotate*pD_old[1];
		pD[1] = sinphi_rotate*pD_old[0] + cosphi_rotate*pD_old[1];
	}
	Grid::normalize_Minkowski<3>(pD,3);
	double kup[4];
	kup[3] = lh->p_kup(com)[3];
	kup[0] = kup[3] * pD[0];
	kup[1] = kup[3] * pD[2];
	kup[2] = kup[3] * pD[1];
	lh->set_p_kup<com>(kup,4);
	//------------------------------------------------------------------------

	// calculate radiation energy in the comoving frame
	const double e_avg = lh->abs_opac(com) > 0 ?
			lh->p_e(com) / (path_length_com * lh->abs_opac(com)) * (1.0 - exp(-lh->abs_opac(com) * path_length_com)) :
			lh->p_e(com);
	Particle fakep;
	LorentzHelper lhtmp(false);
	lhtmp.set_v(lh->velocity(3),3);
	double Dlab[3];

	// deposit amount corresponding to direction actually moved
	fakep = lh->particle_copy(com);
	fakep.e = e_avg;
	fakep.kup[3] = lh->p_kup(com)[3];
	fakep.kup[0] = fakep.kup[3] * Diso[0]; // Diso set above when choosing where to place particle
	fakep.kup[1] = fakep.kup[3] * Diso[1];
	fakep.kup[2] = fakep.kup[3] * Diso[2];
	lhtmp.set_p<com>(&fakep);
	if(randomwalk_n_isotropic <= 0)
		lhtmp.set_distance<com>(path_length_com);
	else
		lhtmp.set_distance<com>(Rcom);
	lhtmp.p_D(lab,Dlab,3);
	zone->distribution[lh->p_s()].count(Dlab, 3, lhtmp.p_nu(com), lhtmp.p_e(lab) * lhtmp.distance(lab));

	// deposit isotropic component
	for(int ip=0; ip<randomwalk_n_isotropic; ip++){
		double Diso_tmp[3];
		grid->isotropic_direction(Diso_tmp,3,&rangen);
		fakep = lh->particle_copy(com);
		fakep.e = e_avg / (double)(randomwalk_n_isotropic);
		fakep.kup[3] = lh->p_kup(com)[3];
		fakep.kup[0] = fakep.kup[3] * Diso_tmp[0];
		fakep.kup[1] = fakep.kup[3] * Diso_tmp[1];
		fakep.kup[2] = fakep.kup[3] * Diso_tmp[2];
		lhtmp.set_p<com>(&fakep);
		lhtmp.set_distance<com>(path_length_com - Rcom);
		lhtmp.p_D(lab,Dlab,3);
		zone->distribution[lh->p_s()].count(Dlab, 3, lhtmp.p_nu(com), lhtmp.p_e(lab) * lhtmp.distance(lab));
	}

	// move the particle to the edge of the sphere
	double xnew[4];
	for(int i=0; i<3; i++) xnew[i] = lh->p_xup()[i] + d3lab[i];
	xnew[3] = lh->p_xup()[3] + lh->distance(lab);
	lh->set_p_xup(xnew,4);
	if(ratio_deposited > 0) lh->scale_p_e( 1.0 - ratio_deposited );
}
