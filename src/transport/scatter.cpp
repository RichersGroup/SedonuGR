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


// decide whether to kill a particle
void Transport::window(EinsteinHelper *eh) const{
	PRINT_ASSERT(eh->p.N,>=,0);
	PRINT_ASSERT(eh->p.fate,!=,rouletted);

	// Roulette if too low energy
	while(eh->p.N<=min_packet_number && eh->p.fate==moving){
		if(rangen.uniform() < 0.5) eh->p.fate = rouletted;
		else eh->p.N *= 2.0;
	}

	if(eh->p.fate == moving){
		PRINT_ASSERT(eh->p.N,>=,min_packet_number);
		PRINT_ASSERT(eh->p.N,<,INFINITY);
		PRINT_ASSERT(eh->p.N,>,0);
	}
}

// choose which type of scattering event to do
void Transport::scatter(EinsteinHelper *eh) const{
	// store the old direction
	double Nold = eh->p.N;
	double kup_tet_old[4];
	for(unsigned i=0; i<4; i++) kup_tet_old[i] = eh->kup_tet[i];
	PRINT_ASSERT(kup_tet_old[3],==,eh->kup_tet[3]);
	
	// sample outgoing energy and set the post-scattered state
	if(use_scattering_kernels) sample_scattering_final_state(eh,kup_tet_old);
	else{
		// sample new direction
		double kup_tet[4];
		isotropic_kup_tet(eh->nu(),kup_tet,&rangen);
		eh->set_kup_tet(kup_tet);
	}

	for(unsigned i=0; i<4; i++){
		grid->fourforce_abs[eh->z_ind][i] += (kup_tet_old[i]*Nold - eh->kup_tet[i]*eh->p.N);
	}
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
			double tmp = 2.0 * exp(-x * (n*pc::pi)*(n*pc::pi)/3.0);
			if(n%2 == 0) tmp *= -1;
			sum += tmp;
	    }

		randomwalk_diffusion_time.set(i,1.0-sum);
	}
	randomwalk_diffusion_time.normalize();
}

//----------------------
// Do a random walk step
//----------------------
void Transport::random_walk(EinsteinHelper *eh) const{
	PRINT_ASSERT(eh->scatopac,>,0);
	PRINT_ASSERT(eh->absopac,>=,0);
	PRINT_ASSERT(eh->p.N,>=,0);
	PRINT_ASSERT(eh->nu(),>=,0);

	const double Rcom = eh->ds_com;
	const double D = eh->scatopac / (3.*pc::c);

	// sample the distance travelled during the random walk
	double path_length_com = pc::c * Rcom*Rcom / D * randomwalk_diffusion_time.invert(rangen.uniform(),&randomwalk_xaxis,-1);
	path_length_com = max(path_length_com,Rcom);

	// move along with the fluid
	double dtau = path_length_com / pc::c;
	for(unsigned i=0; i<4; i++) eh->p.xup[i] += dtau * eh->u[i];

	// determine the average and final neutrino numbers
	double opt_depth = eh->absopac * path_length_com;
	double Nfinal = eh->p.N * exp(-opt_depth);
	double Naverage = (eh->p.N - Nfinal) / (opt_depth);

	// select a random direction
	double kup_tet[4];
	isotropic_kup_tet(eh->nu(),kup_tet,&rangen);
	eh->set_kup_tet(kup_tet);
	eh->ds_com = Rcom;
	eh->p.N = Naverage;
	if(Rcom < INFINITY) tally_radiation(eh);
	move(eh);
	eh->p.N = Nfinal;

	// select a random outward direction
	double kup_tet_final[4] = {0,0,0, kup_tet[3]};
	do{
		isotropic_direction(kup_tet_final,&rangen);
	} while(Metric::dot_Minkowski<3>(kup_tet_final,kup_tet) < 0);
	for(unsigned i=0; i<3; i++) kup_tet_final[i] *= kup_tet_final[3];
	eh->set_kup_tet(kup_tet_final);

	// contribute energy isotropically
	double Eiso;
	Eiso = pc::h*eh->nu() * Naverage * (path_length_com - Rcom);
	grid->distribution[eh->p.s]->add_isotropic(eh->dir_ind, Eiso);
}

//-------------------------------------------------------------
// Sample outgoing neutrino direction and energy
//-------------------------------------------------------------
bool reject_direction(const double mu, const double delta, ThreadRNG* rangen){
	// uniform in direction
	if(delta==0) return false;

	// highly anisotropic - must cut off PDF. Reject if outside bounds
	double min_mu = delta> 1.0 ? delta-2. : -1.0;
	double max_mu = delta<-1.0 ? delta+2. :  1.0;
	if(mu<min_mu || mu>max_mu) return true;

	// If inside bounds, use linear PDF
	double delta_eff = max(min(delta, 1.0), -1.0);
	double mubar = 0.5*(min_mu+max_mu);
	double pdfval = 0.5 + delta_eff * (mu-mubar) / (max_mu-min_mu);
	PRINT_ASSERT(pdfval,<=,1.);
	PRINT_ASSERT(pdfval,>=,0.);
	if(rangen->uniform() > pdfval) return true;
	else return false;
}
void Transport::sample_scattering_final_state(EinsteinHelper *eh, const double kup_tet_old[4]) const{
	PRINT_ASSERT(use_scattering_kernels,>,0);
	PRINT_ASSERT(grid->scattering_delta[eh->p.s].size(),>,0);
	PRINT_ASSERT(grid->scattering_phi0[eh->p.s].size(),>,0);
	PRINT_ASSERT(kup_tet_old[3],==,eh->kup_tet[3]);

	// get spatial component of directional indices
	unsigned dir_ind[NDIMS+2];
	double hyperloc[NDIMS+2];
	for(unsigned i=0; i<NDIMS; i++){
		hyperloc[i] = eh->p.xup[i];
		dir_ind[i] = eh->dir_ind[i];
	}
	dir_ind[NDIMS] = eh->dir_ind[NDIMS];

	// get outgoing energy bin w/ rejection sampling
	unsigned igout;
	double P, phi0avg, nubar;
	do{
		igout = rangen.uniform_discrete(0, grid->nu_grid_axis.size()-1);
		dir_ind[NDIMS+1] = igout;
		nubar = 0.5 * (grid->nu_grid_axis.top[igout] + grid->nu_grid_axis.bottom(igout));
		hyperloc[NDIMS] = grid->nu_grid_axis.mid[dir_ind[NDIMS]];
		hyperloc[NDIMS+1] = nubar;
		phi0avg = grid->scattering_phi0[eh->p.s].interpolate(hyperloc,dir_ind);
		P = phi0avg * grid->nu_grid_axis.delta(igout) / grid->scat_opac[eh->p.s][eh->eas_ind];
		PRINT_ASSERT(P,<=,1.0);
	} while(rangen.uniform() > P);

	// uniformly sample within zone
	unsigned global_index = grid->scattering_phi0[eh->p.s].direct_index(dir_ind);
	double out_nu = rangen.uniform(grid->nu_grid_axis.bottom(igout), grid->nu_grid_axis.top[igout]);
	double phi_interpolated = phi0avg; // + grid->scattering_phi0[eh->p.s].dydx[global_index][NDIMS+1][0]*(out_nu - nubar);
	eh->p.N *= phi_interpolated / phi0avg;

	// get scattering delta
	hyperloc[NDIMS] = eh->nu();
	hyperloc[NDIMS+1] = out_nu;
	double delta = grid->scattering_delta[eh->p.s].interpolate(hyperloc,dir_ind);
	PRINT_ASSERT(fabs(delta),<,3.0);

	// sample the new direction, but only if not absurdly forward/backward peaked
	// (delta=2.8 corresponds to a possible factor of 10 in the neutrino weight)
	if(fabs(delta) < 2.8){
		double kup_tet_new[4];
		double costheta=0;
		do{
			isotropic_kup_tet(out_nu, kup_tet_new, &rangen);
			costheta = Metric::dot_Minkowski<3>(kup_tet_new,kup_tet_old) / (kup_tet_old[3]*kup_tet_new[3]);
		} while(reject_direction(costheta, delta, &rangen));
		eh->set_kup_tet(kup_tet_new);
		PRINT_ASSERT(eh->p.N,<,1e99);
	}
	else{
		eh->scale_p_frequency(out_nu/eh->nu());
	}
}
