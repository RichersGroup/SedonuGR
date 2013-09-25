#pragma warning disable 161
#include <limits>
#include <omp.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include "transport.h"
#include "physical_constants.h"
#include "species_general.h"

namespace pc = physical_constants;


void transport::propagate_particles(double dt)
{
  vector<int> n_active(species_list.size(),0);
  vector<int> n_escape(species_list.size(),0);
  double e_esc = 0;
  double N;
  #pragma omp parallel shared(n_active,n_escape,e_esc,N) firstprivate(dt) 
  {

    //--- MOVE THE PARTICLES AROUND ---
    #pragma omp for schedule(guided) reduction(+:e_esc)
    for(int i=0; i<particles.size(); i++){
      particle* p = &particles[i];
      #pragma omp atomic
      n_active[p->s]++;
      propagate(p,dt);
      if(p->fate == escaped){
	#pragma omp atomic
	n_escape[p->s]++;
	species_list[p->s]->spectrum.count(p->t, p->nu, p->e, p->D);
	e_esc += p->e;
      }
    } //implied barrier

    //--- REMOVE THE DEAD PARTICLES ---
    #pragma omp single nowait
    {
      vector<particle>::iterator pIter = particles.begin();
      while(pIter != particles.end()){
	if(pIter->fate==absorbed || pIter->fate==escaped){
	  *pIter = particles[particles.size()-1];
	  particles.pop_back();
	}
	else pIter++;
      }
    }

    //--- DETERMINE THE NORMALIZATION FACTOR ---
    // accounts for particles scattering back to core / being absorbed
    // and the fact that emitted particles are given a specific energy in their
    // rest frame, which makes the lab frame emitted energy different from L_net.
    #pragma omp single
    {
      N = 1;
      if(iterate){
	if(e_esc>0) N = L_net / e_esc;
	else if(verbose) cout << "# WARNING: no energy escaped. Setting normalization to 1." << endl;
      }
    }

    //--- NORMALIZE THE GRID QUANTITIES ---
    #pragma omp for nowait
    for(int i=0; i<grid->z.size(); i++){
      grid->z[i].e_rad *= N;
      grid->z[i].e_abs *= N;
      grid->z[i].l_abs *= N;
    }
  } //#pragma omp parallel
  
  //--- OUPUT ESCAPE STATISTICS AND NORMALIZE SPECTRUM ---
  for(int i=0; i<species_list.size(); i++){
    if(n_escape[i]>0) species_list[i]->spectrum.rescale(N);
    double per_esc = (100.0*n_escape[i])/n_active[i];
    if (verbose && iterate){
      if(n_active[i]>0) printf("# %i/%i %s escaped. (%3.2f%%)\n", n_escape[i], n_active[i], species_list[i]->name.c_str(), per_esc);
      else printf("# No active %s.\n", species_list[i]->name.c_str());
    }
  }
}

//--------------------------------------------------------
// Propagate a single monte carlo particle until
// it  escapes, is absorbed, or the time step ends
//--------------------------------------------------------
void transport::propagate(particle* p, double dt)
{
  enum ParticleEvent {scatter, boundary, tstep};
  ParticleEvent event;

  if((p->ind < 0) || (p->ind > grid->z.size())){
    cout << "ERROR: particle coming in to propagate has invalid index" << endl;
    exit(3);
  }

  p->fate = moving;

  // time of end of timestep
  double tstop = t_now + dt;

  // local variables
  zone* zone;                                      // pointer to current zone
  double tau_r;                                    // random optical depth
  double d_bn,d_sc,d_tm,this_d;                    // interaction distances
  double dshift;                                   // doppler shift
  double opac, abs_frac, opac_lab;                 // opacity variables
  double this_E, this_E_comoving, this_l_comoving; // for calculating radiation energy and energy/lepton number absorbed
  double z,z2;                                     // random numbers

  // propagate until this flag is set
  while (p->fate == moving)
  {
    // set pointer to current zone
    zone = &(grid->z[p->ind]);

    // maximum step size inside zone
    d_bn = step_size * grid->zone_min_length(p->ind);

    // doppler shift from comoving to lab
    dshift = dshift_comoving_to_lab(p);

    // get local opacity and absorption fraction
    species_list[p->s]->get_opacity(p,dshift,&opac,&abs_frac);

    // convert opacity from comoving to lab frame for the purposes of
    // determining the interaction distance in the lab frame
    // This corresponds to equation 90.8 in Mihalas&Mihalas. You multiply
    // the comoving opacity by nu_0 over nu, which is why you
    // multiply by dshift instead of dividing by dshift here
    opac_lab = opac*dshift;

    // random optical depth to next interaction
    tau_r = -1.0*log(1 - rangen.uniform());

    // step size to next interaction event
    d_sc  = tau_r/opac_lab;
    if (opac_lab == 0) d_sc = numeric_limits<double>::infinity();
    if (d_sc < 0){
      cout << "ERROR: negative interaction distance!\n" << endl;
      cout << __FILE__ << ":" << __LINE__ << endl;
      exit(15);
    }

    // find distance to end of time step
    d_tm = (tstop - p->t)*pc::c;
    if (iterate) d_tm = numeric_limits<double>::infinity(); // i.e. let all particles escape

    // find out what event happens (shortest distance)
    if ( (d_sc < d_bn) && (d_sc < d_tm) ){
      event  = scatter;
      this_d = d_sc;
    }
    else if (d_bn < d_tm){
      event  = boundary;
      this_d = d_bn;
    }
    else{
      event  = tstep;
      this_d = d_tm;
    }

    // tally in contribution to zone's radiation energy (both *lab* frame)
    this_E = p->e*this_d;
    #pragma omp atomic
    zone->e_rad += this_E;

    // store absorbed energy in *comoving* frame
    // (will turn into rate by dividing by dt later)
    // Extra dshift definitely needed here (two total)
    // to convert both p->e and this_d to the comoving frame
    this_E_comoving = this_E * dshift * dshift;
    #pragma omp atomic
    zone->e_abs += this_E_comoving * (opac*abs_frac*zone->eps_imc);

    // store absorbed lepton number (same in both frames, except for the
    // factor of this_d which is divided out later
    if(species_list[p->s]->lepton_number != 0){
      this_l_comoving = species_list[p->s]->lepton_number * p->e/(p->nu*pc::h) * this_d*dshift;
      #pragma omp atomic
      zone->l_abs += this_l_comoving * (opac*abs_frac*zone->eps_imc);
    }

    // TODO - put back in radiation force tally here
    // fx_rad =

    // move particle the distance
    p->x[0] += this_d*p->D[0];
    p->x[1] += this_d*p->D[1];
    p->x[2] += this_d*p->D[2];

    // advance the time
    p->t = p->t + this_d/pc::c;

    // ---------------------------------
    // Do if scatter
    // ---------------------------------
    if (event == scatter)
    {
      // random number to check for scattering or absorption
      z = rangen.uniform();

      // decide whether to scatter or absorb
      if (z > abs_frac) isotropic_scatter(p,0);
      else
      {
	// RNG for deciding whether to do effective scattering
	if (radiative_eq) z2 = 2;
	else z2 = rangen.uniform();

	// choose whether to do an effective scatter (i.e. particle is absorbed
	// but, since we require energy in = energy out it is re-emitted)
	// otherwise, really absorb (kill) it
	if (z2 > zone->eps_imc) isotropic_scatter(p,1);
	else p->fate = absorbed;
      }
    }

    // ---------------------------------
    // do if time step end
    // ---------------------------------
    else if (event == tstep) p->fate = stopped;

    // Find position of the particle now
    p->ind = grid->get_zone(p->x);
    if (p->ind == -1) p->fate = absorbed;
    if (p->ind == -2) p->fate = escaped;

    // check for inner boundary absorption
    if (p->r() < r_core) p->fate = absorbed;
  }
}
