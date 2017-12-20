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
void Transport::event_interact(LorentzHelper* lh, int *z_ind){
	PRINT_ASSERT(*z_ind,>=,0);
	PRINT_ASSERT(*z_ind,<,(int)grid->z.size());
	PRINT_ASSERT(lh->abs_fraction(),>=,0.0);
	PRINT_ASSERT(lh->abs_fraction(),<=,1.0);
	PRINT_ASSERT(lh->p_N(),>,0);
	PRINT_ASSERT(lh->p_fate(),==,moving);

	// absorb part of the packet
	if(!exponential_decay) lh->scale_p_number(1.0 - lh->abs_fraction());
	scatter(lh, z_ind);

	// resample the path length
	if(lh->p_fate()==moving) lh->set_p_tau(sample_tau(&rangen));

	// window the particle
	if(lh->p_fate()==moving) window(lh,*z_ind);

	// sanity checks
	if(lh->p_fate()==moving){
		PRINT_ASSERT(lh->p_nu(),>,0);
		PRINT_ASSERT(lh->p_nu(),>,0);
	}
}


// decide whether to kill a particle
void Transport::window(LorentzHelper *lh, const int z_ind){
	PRINT_ASSERT(lh->p_N(),>=,0);
	PRINT_ASSERT(lh->p_fate(),!=,rouletted);

	// Roulette if too low energy
	while(lh->p_N()<=min_packet_number && lh->p_fate()==moving){
		if(rangen.uniform() < 0.5) lh->set_p_fate(rouletted);
		else lh->scale_p_number(2.0);
	}
	if(lh->p_fate()==moving) PRINT_ASSERT(lh->p_N(),>=,min_packet_number);

	// split if too high energy, if enough space, and if in important region
	double ratio = lh->p_N() / max_packet_number;
	int n_new = (int)ratio;
	if(ratio>1.0 && (int)particles.size()+n_new<max_particles){
		lh->scale_p_number( 1.0 / (double)(n_new+1) );
		for(int i=0; i<n_new; i++){
			#pragma omp critical
			particles.push_back(lh->particle_copy(lab));
		}
	}

	if(lh->p_fate() == moving){
		PRINT_ASSERT(lh->p_N(),<,INFINITY);
		PRINT_ASSERT(lh->p_N(),>,0);
	}
	if((int)particles.size()>=max_particles && verbose && rank0){
		cout << "max_particles: " << max_particles << endl;
		cout << "particles.size(): " << particles.size() << endl;
		cout << "WARNING: max_particles is too small to allow splitting." << endl;
	}
}

// choose which type of scattering event to do
void Transport::scatter(LorentzHelper *lh, int *z_ind) const{
	PRINT_ASSERT(lh->abs_fraction(),>=,0);
	PRINT_ASSERT(lh->abs_fraction(),<=,1.0);
	PRINT_ASSERT(lh->net_opac(com),>=,0);
	bool did_random_walk = false;

	// try to do random walk approximation in scattering-dominated diffusion regime
	if(randomwalk_sphere_size>0){

		double D = pc::c / (3.0 * lh->scat_opac(com)); // diffusion coefficient (cm^2/s)

		// if the optical depth is below our threshold, don't do random walk
		// (first pass to avoid doing lots of math)
		double Rlab_min = randomwalk_sphere_size * grid->zone_min_length(*z_ind);
		double Rlab_boundary = grid->zone_cell_dist(lh->p_xup(),*z_ind);
		double Rlab = max(Rlab_min, Rlab_boundary);
		if(lh->scat_opac(com) * Rlab >= randomwalk_min_optical_depth){
			// determine maximum comoving sphere size
			const double* v = lh->velocity(3);
			double vabs = sqrt(Grid::dot_Minkowski<3>(v,v));
			double gamma = LorentzHelper::lorentz_factor(v,3);

			double Rcom = 0;
			if(Rlab==0) Rcom = 0;
			else if(Rlab==INFINITY) Rcom = randomwalk_sphere_size * randomwalk_min_optical_depth / (lh->abs_opac(com)>0 ? lh->abs_opac(com) : lh->scat_opac(com));
			else Rcom =  2. * Rlab / gamma / (1. + sqrt(1. + 4.*Rlab*vabs*randomwalk_max_x / (gamma*D) ) );

			// if the optical depth is below our threshold, don't do random walk
			if(lh->scat_opac(com) * Rcom >= randomwalk_min_optical_depth){
				random_walk(lh, Rcom, D, *z_ind);
				*z_ind = grid->zone_index(lh->p_xup());\
				boundary_conditions(lh, z_ind);
				//if(new_ind==-2) lh->set_p_fate(escaped);
				did_random_walk = true;
			}
		}
	}

	// isotropic scatter if can't do random walk
	if(!did_random_walk && lh->p_fate()==moving){
		// store the old direction
		double kup_old[3];
		kup_old[0] = lh->p_kup(com)[0];
		kup_old[1] = lh->p_kup(com)[1];
		kup_old[2] = lh->p_kup(com)[2];

		// sample new direction
		double kup[4];
		grid->isotropic_kup_tet(lh->p_nu(),kup,lh->p_xup(),&rangen);
		lh->set_p_kup<com>(kup,4);

		// get the dot product between the old and new directions
		double cosTheta = grid->dot3(kup,kup_old,lh->p_xup()) / (lh->p_nu() * lh->p_nu() * 4. * pc::pi * pc::pi / (pc::c * pc::c));
		PRINT_ASSERT(fabs(cosTheta),<=,1.0);

		// sample outgoing energy and set the post-scattered state
		if(use_scattering_kernels){
			const double tmp_nu = lh->p_nu();
			sample_scattering_final_state(*z_ind,*lh,cosTheta);
			double dep_energy = (tmp_nu - lh->p_nu()) * lh->p_N()*pc::h;
			if(dep_energy>0) grid->z[*z_ind].e_abs  += dep_energy;
			else             grid->z[*z_ind].e_emit -= dep_energy;
		}
	}
}



//---------------------------------------------------------------------
// Randomly select an optical depth through which a particle will move.
//---------------------------------------------------------------------
double Transport::sample_tau(ThreadRNG* rangen){
	double tau;

	do{
		tau = -log(rangen->uniform());
	} while(tau >= INFINITY);

	return tau;
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
void Transport::random_walk(LorentzHelper *lh, const double Rcom, const double D, const int z_ind) const{
	PRINT_ASSERT(lh->scat_opac(com),>,0);
	PRINT_ASSERT(lh->abs_opac(com),>=,0);
	PRINT_ASSERT(lh->p_N(),>=,0);
	PRINT_ASSERT(lh->p_nu(),>=,0);

	// set pointer to the current zone
	Zone* zone;
	zone = &(grid->z[z_ind]);

	// sample the distance travelled during the random walk
	double path_length_com = pc::c * Rcom*Rcom / D * randomwalk_diffusion_time.invert(rangen.uniform(),&randomwalk_xaxis,-1);
	//PRINT_ASSERT(path_length_com,>=,Rcom);
	path_length_com = max(path_length_com,Rcom);

	//=================================//
	// pick a random direction to move //
	//=================================//
	double Diso[3], xnew[4];
	grid->isotropic_direction(Diso,&rangen);
	double displacement4[4] = {Rcom*Diso[0], Rcom*Diso[1], Rcom*Diso[2], path_length_com};
	double d3com[3] = {displacement4[0], displacement4[1], displacement4[2]};
	LorentzHelper::transform_cartesian_4vector_c2l(zone->u, displacement4);
	double d3lab[3] = {displacement4[0], displacement4[1], displacement4[2]};
	xnew[0] = lh->p_xup()[0] + d3lab[0];
	xnew[1] = lh->p_xup()[1] + d3lab[1];
	xnew[2] = lh->p_xup()[2] + d3lab[2];
	xnew[3] = lh->p_xup()[3] + lh->distance(lab);

	//===============================//
	// pick a random final direction //
	//===============================//
	double pD[3];
	double outward_phi = 2*pc::pi * rangen.uniform();
	double outward_mu = rangen.uniform();
	pD[0] = outward_mu*cos(outward_phi);
	pD[1] = outward_mu*sin(outward_phi);
	pD[2] = sqrt(1.0 - outward_mu*outward_mu);

	// get the displacement vector polar coordinates
	// theta_rotate is angle away from z-axis, not from xy-plane
	Grid::normalize_Minkowski<3>(d3com);
	double costheta_rotate = d3com[2];
	double sintheta_rotate = sqrt(1.0 - d3com[2]*d3com[2]);

	if(abs(sintheta_rotate) < TINY) pD[2] *= costheta_rotate>0 ? 1.0 : -1.0;
	else{

		// first rotate away from the z axis along y=0 (move it toward x=0)
		Grid::normalize_Minkowski<3>(pD);
		double pD_old[3];
		for(int i=0; i<3; i++) pD_old[i] = pD[i];
		pD[0] =  costheta_rotate*pD_old[0] + sintheta_rotate*pD_old[2];
		pD[2] = -sintheta_rotate*pD_old[0] + costheta_rotate*pD_old[2];

		// second rotate around the z axis, away from the x-axis
		double cosphi_rotate = d3com[0] / sintheta_rotate;
		double sinphi_rotate = d3com[1] / sintheta_rotate;
		Grid::normalize_Minkowski<3>(pD);
		for(int i=0; i<3; i++) pD_old[i] = pD[i];
		pD[0] = cosphi_rotate*pD_old[0] - sinphi_rotate*pD_old[1];
		pD[1] = sinphi_rotate*pD_old[0] + cosphi_rotate*pD_old[1];
	}
	Grid::normalize_Minkowski<3>(pD);
	double kup_new[4];
	kup_new[3] = lh->p_kup(com)[3];
	kup_new[0] = kup_new[3] * pD[0];
	kup_new[1] = kup_new[3] * pD[2];
	kup_new[2] = kup_new[3] * pD[1];

	//==============================//
	// setup for radiation tallying //
	//==============================//
	LorentzHelper lhtmp = *lh;
	lhtmp.exponential_decay = true;
	lhtmp.set_distance<com>(path_length_com);
	Particle fakeP = lh->particle_copy(com);

	//===============================//
	// DEPOSIT DIRECTIONAL COMPONENT //
	//===============================//

	fakeP.N = lh->p_N();
	if(randomwalk_n_isotropic > 0)
		fakeP.N *= Rcom / path_length_com;
	PRINT_ASSERT(fakeP.N,>=,0);
	fakeP.kup[0] = fakeP.kup[3] * Diso[0]; // Diso set above when choosing where to place particle
	fakeP.kup[1] = fakeP.kup[3] * Diso[1];
	fakeP.kup[2] = fakeP.kup[3] * Diso[2];
	lhtmp.set_p<com>(&fakeP);
	tally_radiation(&lhtmp,z_ind);

	//=============================//
	// DEPOSIT ISOTROPIC COMPONENT //
	//=============================//

	if(randomwalk_n_isotropic > 0){
		fakeP.N = lh->p_N() * (path_length_com - Rcom)/path_length_com / (double)randomwalk_n_isotropic;
		if(fakeP.N > 0) for(int ip=0; ip<randomwalk_n_isotropic; ip++){
			// select a random direction
			double Diso_tmp[3];
			grid->isotropic_direction(Diso_tmp,&rangen);
			fakeP.kup[0] = fakeP.kup[3] * Diso_tmp[0];
			fakeP.kup[1] = fakeP.kup[3] * Diso_tmp[1];
			fakeP.kup[2] = fakeP.kup[3] * Diso_tmp[2];

			// tally the contribution
			lhtmp.set_p<com>(&fakeP);
			tally_radiation(&lhtmp,z_ind);
		}
	}

	//=============================================//
	// move the particle to the edge of the sphere //
	//=============================================//
	lh->set_p_xup(xnew,4);
	lh->set_p_kup<com>(kup_new,4);
	lh->scale_p_number( exp(-lh->abs_opac(com) * path_length_com) );
}

//-------------------------------------------------------------
// Sample outgoing neutrino direction and energy
//-------------------------------------------------------------
void Transport::sample_scattering_final_state(const int z_ind, LorentzHelper &lh, const double cosTheta) const{
	PRINT_ASSERT(use_scattering_kernels,>,0);
	PRINT_ASSERT(species_list[lh.p_s()]->scattering_delta.size(),>,0);
	PRINT_ASSERT(species_list[lh.p_s()]->normalized_phi0.size(),>,0);

	// get ingoing energy index
	int igin = grid->nu_grid.locate(lh.p_nu());
	if(igin == grid->nu_grid.size()) igin--;

	// get outgoing energy
	double out_nu = species_list[lh.p_s()]->normalized_phi0[z_ind][igin].invert(rangen.uniform(),&grid->nu_grid);
	lh.scale_p_energy(out_nu/lh.p_nu());

	// get outgoing energy index
	int igout = grid->nu_grid.locate(out_nu);

	// bias outgoing direction to be isotropic. Very inefficient for large values of delta.
	double delta = species_list[lh.p_s()]->scattering_delta[z_ind][igin][igout];
	PRINT_ASSERT(fabs(delta),<,3.0);
	if(fabs(delta)<=1.0) lh.scale_p_number(1.0 + delta*cosTheta);
	else{
		double b = 2.*fabs(delta) / (3.-fabs(delta));
		if(delta>1.0) lh.scale_p_number( pow(1.+cosTheta, b) );
		else          lh.scale_p_number( pow(1.-cosTheta, b) );
	}


	// sample outgoing direction
	//double U = sim->rangen.uniform();
	//double tiny = 1.e-4;
	//if(abs(delta)<1.e-4) // isotropic
	//	*outCosTheta = 2.*U - 1.;
	//else{
	//	double b2m4ac = sqrt(1. - delta*(2.-delta-4.*U));
	//	double sgn = (delta<0 ? -1.0 : 1.0);
	//	*outCosTheta = 1./delta * (-1. + sgn*b2m4ac);
	//}
}
