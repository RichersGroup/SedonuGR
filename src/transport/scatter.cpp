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
	PRINT_ASSERT(eh->N,>=,0);
	PRINT_ASSERT(eh->fate,!=,rouletted);

	// Roulette if too low energy
	if(eh->N == 0) eh->fate = rouletted;
	else while(eh->N/eh->N0 < min_packet_weight and eh->fate==moving){
		if(rangen.uniform() < 0.5) eh->fate = rouletted;
		else eh->N *= 2.0;
	}

	if(eh->fate == moving){
		PRINT_ASSERT(eh->N/eh->N0,>=,min_packet_weight);
		PRINT_ASSERT(eh->N,<,INFINITY);
		PRINT_ASSERT(eh->N,>,0);
	}
}

// choose which type of scattering event to do
void Transport::scatter(EinsteinHelper *eh, const ParticleEvent event) const{
	assert(event==elastic_scatter or event==inelastic_scatter);

	// store the old direction
	double Nold = eh->N;
	Tuple<double,4> kup_tet_old = eh->kup_tet;
	
	// sample outgoing energy and set the post-scattered state
	if(event==inelastic_scatter){
		PRINT_ASSERT(eh->inelastic_scatopac,>=,0);
		if(eh->inelastic_scatopac > 0)
			sample_scattering_final_state(eh,kup_tet_old);
	}
	else if(event==elastic_scatter){
		PRINT_ASSERT(eh->scatopac,>=,0);
		if(eh->scatopac > 0){
			Tuple<double,4> kup_tet = eh->kup_tet;
			isotropic_kup_tet(kup_tet,&rangen);
			eh->set_kup_tet(kup_tet);
		}
	}

	grid->fourforce_abs[eh->z_ind] += (kup_tet_old - eh->kup_tet) * eh->N / eh->zone_fourvolume;
}

double Pescape(double x, int sumN){
  double sum = 0;
  for(int n=1; n<=sumN; n++){
    double tmp = 2.0 * exp(-x * (n*pc::pi)*(n*pc::pi));
    if(n%2 == 0) tmp *= -1;
    sum += tmp;
  }
  return 1.0-sum;
}


//-------------------------------------------------------
// Initialize the CDF that determines particle dwell time
// result is D*t/(R^2)
//-------------------------------------------------------
void Transport::init_randomwalk_cdf(Lua* lua){
	randomwalk_sumN = lua->scalar<int>("randomwalk_sumN");
	int npoints = lua->scalar<int>("randomwalk_npoints");
	randomwalk_max_x = lua->scalar<double>("randomwalk_max_x");

	randomwalk_diffusion_time.resize(npoints);
	randomwalk_diffusion_time.interpolation_order = 1;
	randomwalk_xaxis = Axis(0,randomwalk_max_x,npoints);

	#pragma omp parallel for
	for(int i=1; i<=npoints; i++)
	  randomwalk_diffusion_time.set(i,Pescape(randomwalk_xaxis.top[i], randomwalk_sumN));
	randomwalk_diffusion_time.normalize();
}

//----------------------
// Do a random walk step
//----------------------
void Transport::random_walk(EinsteinHelper *eh) const{
	PRINT_ASSERT(eh->scatopac,>,0);
	PRINT_ASSERT(eh->absopac,>=,0);
	PRINT_ASSERT(eh->N,>=,0);
	PRINT_ASSERT(eh->nu(),>=,0);

	// invert the randomwalk CDF
	const double Rcom = eh->ds_com;
	const double D = pc::c / (3.*eh->scatopac);
	double chi_esc = D/(pc::c*Rcom); // using t=Rcom/c
	double Pesc_cdf  = randomwalk_diffusion_time.interpolate_cdf(chi_esc,&randomwalk_xaxis);
	double Pesc_real = -exp(-eh->scatopac*Rcom);
	double U = rangen.uniform();
	double path_length_com = 0;
	if(U<Pesc_real) path_length_com = Rcom;
	else{
		// make U go from Pesc_cdf to 1 (following stretching in Foucart2018)
		double a = (1. - Pesc_cdf) / (1. - Pesc_real);
		U = (U-Pesc_real) * a + Pesc_cdf;
		PRINT_ASSERT(U,>=,0-TINY);
		PRINT_ASSERT(U,<=,1+TINY);
		U = max(min(U,1.),0.);

		// invert the randomwalk CDF
		double chi = randomwalk_diffusion_time.invert(U,&randomwalk_xaxis,-1);
		path_length_com = pc::c * Rcom*Rcom / D * chi;
		PRINT_ASSERT(path_length_com,>=,Rcom*(1.-TINY));
		path_length_com = max(path_length_com,Rcom);
	}

	// foucart 2018 description
	double f_free = Rcom/path_length_com;
	PRINT_ASSERT(f_free,>=,0);
	PRINT_ASSERT(f_free,<=,1);
	double ds_free = path_length_com * f_free;
	double ds_fl   = path_length_com * (1.-f_free);

	// Foucart2018 31-33
	double ds_adv = 0;
	double A=0,B=0;
	if(Metric::dot_Minkowski<3>(eh->u,eh->u) > 1e-10){
	  double gtt = (DO_GR ? eh->g.get(3,3) : -1. );
	  Tuple<double,4> ulow = eh->g.lower<4>(eh->u);
	  double A_B = -(ulow[3] + sqrt(ulow[3]*ulow[3] + gtt)) / gtt;
	  B = eh->kup_tet[3] / (1. - A_B*ulow[3]);
	  A = A_B * B;
	  double up_up = (A+B*eh->u[3]) / (B*eh->u[3]);
	  ds_adv = ds_fl * up_up;
	  PRINT_ASSERT(ds_adv,<=,ds_fl);
	  PRINT_ASSERT(ds_adv,>=,0);
	}
	double ds_iso = ds_fl - ds_adv;

	//================//
	// Isotropic Step //
	//================//
	if(ds_iso>0){
	  // determine the average and final neutrino numbers
	  double Naverage = eh->N, Nfinal = eh->N, Nold = eh->N;
	  if(eh->absopac > 0){
	    double opt_depth = eh->absopac * ds_iso;
	    Nfinal = eh->N * exp(-opt_depth);
	    if((eh->N-Nfinal)/eh->N < TINY)
	      Naverage = (eh->N + Nfinal) / 2.;
	    else
	      Naverage = (eh->N - Nfinal) / (opt_depth);
	  }
	  PRINT_ASSERT(Naverage,<=,Nold);
	  PRINT_ASSERT(Nfinal,<=,Naverage);
	  
	  // contribute isotropically
	  double Eiso = eh->kup_tet[3] * Naverage * ds_iso / (eh->zone_fourvolume*pc::c);
	  grid->distribution[eh->s]->add_isotropic_single(eh->dir_ind, Eiso);
	  grid->l_abs[eh->z_ind] += (Nold - Nfinal) * species_list[eh->s]->lepton_number / eh->zone_fourvolume;
	  grid->fourforce_abs[eh->z_ind] += eh->kup_tet * (Nold - Nfinal) / eh->zone_fourvolume;
	  
	  // move neutrino forward in time
	  eh->xup[3] += ds_iso * eh->u[3];
	  eh->N = Nfinal;
	  window(eh);
	}
	  
	//================//
	// Advection step // 
	//================//
	if(ds_adv>0 and eh->fate==moving){
	  // set the momentum in the direction of the fluid
	  Tuple<double,4> tup;
	  tup[0] = tup[1] = tup[2] = 0;
	  tup[3] = 1;
	  Tuple<double,4> kup_tet_old = eh->kup_tet;
	  eh->kup = tup*A + eh->u*B;
	  PRINT_ASSERT(abs(eh->g.dot<4>(eh->kup,eh->kup)),<,TINY);
	  eh->renormalize_kup();
	  PRINT_ASSERT(abs(kup_tet_old[3]-eh->kup_tet[3])/kup_tet_old[3],<,TINY);

	  // account for change in the fluid
	  grid->fourforce_abs[eh->z_ind] += (kup_tet_old - eh->kup_tet) * eh->N / eh->zone_fourvolume;
	  
	  // move for the small timestep
	  eh->ds_com = ds_adv;
	  move(eh);
	}

	//=====================//
	// Free-streaming step //
	//=====================//
	if(ds_free>0 and eh->fate==moving){
	  // select a random direction
	  Tuple<double,4> kup_tet_old = eh->kup_tet;
	  Tuple<double,4> kup_tet = eh->kup_tet;
	  isotropic_kup_tet(kup_tet,&rangen);
	  eh->set_kup_tet(kup_tet);

	  // account for change in the fluid
	  grid->fourforce_abs[eh->z_ind] += (kup_tet_old - eh->kup_tet) * eh->N / eh->zone_fourvolume;

	  // move forward
	  eh->ds_com = ds_free;
	  move(eh);
	  if(eh->fate!=moving) return;

	  // select a random outward direction. Use delta=2 to make pdf=costheta
	  kup_tet_old = eh->kup_tet;
	  kup_tet     = eh->kup_tet;
	  do{
	    isotropic_kup_tet(kup_tet,&rangen);
	  } while(reject_direction(Metric::dot_Minkowski<3>(kup_tet_old,kup_tet)/(kup_tet[3]*kup_tet[3]), 2.) );
	  eh->set_kup_tet(kup_tet);
	
	  // account for change in the fluid
	  grid->fourforce_abs[eh->z_ind] += (kup_tet_old - eh->kup_tet) * eh->N / eh->zone_fourvolume;
	}
}

void Transport::sample_scattering_final_state(EinsteinHelper *eh, const Tuple<double,4>& kup_tet_old) const{
	PRINT_ASSERT(eh->inelastic_scatopac,>,0);
	PRINT_ASSERT(grid->scattering_delta[eh->s].size(),>,0);
	PRINT_ASSERT(kup_tet_old[3],==,eh->kup_tet[3]);

	// rejection sampling to get outgoing frequency bin.
	size_t igout;
	double P;
	do{
		igout = rangen.uniform_discrete(0, grid->nu_grid_axis.size()-1);
		P = grid->partial_scat_opac[eh->s][igout].interpolate(eh->icube_spec) / eh->inelastic_scatopac;
		PRINT_ASSERT(P-1.0,<=,TINY);
		PRINT_ASSERT(P,>=,0.0);
	} while(rangen.uniform() > P);

	// Scatter to the center of the new bin.
	double outnu = grid->nu_grid_axis.mid[igout];
	//cout<<"inelastic"<<endl;
	// interpolate the kernel anisotropy
	double delta = grid->scattering_delta[eh->s][igout].interpolate(eh->icube_spec);
	PRINT_ASSERT(fabs(delta),<,3.0);

	// rejection sample the new direction, but only if not absurdly forward/backward peaked
	// (delta=2.8 corresponds to a possible factor of 10 in the neutrino weight)
	Tuple<double,4> kup_tet_new;
	kup_tet_new[3] = outnu * pc::h;
	if(fabs(delta) < 2.8){
		double costheta=0;
		do{
			isotropic_kup_tet(kup_tet_new, &rangen);
			costheta = Metric::dot_Minkowski<3>(kup_tet_new,kup_tet_old) / (kup_tet_old[3]*kup_tet_new[3]);
		} while(reject_direction(costheta, delta));
	}
	else{
	        kup_tet_new = eh->kup_tet * outnu / eh->nu();
		if(delta<0) for(size_t i=0; i<3; i++) kup_tet_new[i] = -eh->kup_tet[i];
	}

	//check whether scattering should be blocked
	EinsteinHelper eh_new = *eh;
	eh_new.set_kup_tet(kup_tet_new);
	update_eh_k_opac(&eh_new);
	double blocking = grid->fblock[eh->s].interpolate(eh_new.icube_spec);
	if(rangen.uniform() < blocking) return;

	*eh = eh_new;
	PRINT_ASSERT(eh->N,<,1e99);
}

