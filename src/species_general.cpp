#include "species_general.h"
#include <gsl/gsl_rng.h>
#include <math.h>
#include "transport.h"
#include "physical_constants.h"
#include "Lua.h"

namespace pc = physical_constants;

void species_general::init(Lua* lua, transport* simulation)
{
  // set the pointer to see the simulation info
  sim = simulation;

  // call child's init function
  myInit(lua);
}


double species_general::sample_core_nu()
{
  // sample to find the frequency bin to use
  double z = gsl_rng_uniform(sim->rangen);
  int ilam = core_emis.sample(z);

  // sample uniformily in selected frequency bin 
  z = gsl_rng_uniform(sim->rangen);
  return nu_grid.sample(ilam,z);
}

//------------------------------------------------------------
// get the doppler shift from lab to comoving
//------------------------------------------------------------
double species_general::dshift_lab_to_comoving(particle p)
{
  double v[3];
  sim->grid->velocity_vector(p.ind,p.x,v);
  
  // outgoing velocity vector
  double beta = 0.0;
  for (int i=0;i<3;i++)  beta += v[i]*v[i]; 
  beta = sqrt(beta)/pc::c;
  double gamma  = 1.0/sqrt(1 - beta*beta); 

  // doppler shifts outgoing
  double vdp    = (p.D[0]*v[0] + p.D[1]*v[1] + p.D[2]*v[2]);
  double dshift = gamma*(1 - vdp/pc::c);
  return dshift;
}


//------------------------------------------------------------
// get the doppler shift from comoving to lab 
//------------------------------------------------------------
double species_general::dshift_comoving_to_lab(particle p)
{
  double v[3];
  sim->grid->velocity_vector(p.ind,p.x,v);
  
  // outgoing velocity vector
  double beta = 0.0;
  for (int i=0;i<3;i++) 
  {
    v[i] = -1*v[i];
    beta += v[i]*v[i]; 
  }
  beta = sqrt(beta)/pc::c;
  double gamma  = 1.0/sqrt(1 - beta*beta); 

  // doppler shifts outgoing
  double vdp    = (p.D[0]*v[0] + p.D[1]*v[1] + p.D[2]*v[2]);
  double dshift = gamma*(1 - vdp/pc::c);
  return dshift;
}

//------------------------------------------------------------
// do a lorentz transformation; modifies the energy, frequency 
// and direction vector of the particle
// sign =  1 for comoving to lab
// sign = -1 for lab to comoving
//------------------------------------------------------------

void species_general::transform_comoving_to_lab(particle &p)
{
  lorentz_transform(p,1);
}
void species_general::transform_lab_to_comoving(particle &p)
{
  lorentz_transform(p,-1);
}

void species_general::lorentz_transform(particle &p, double sign)
{
  double v[3];
  sim->grid->velocity_vector(p.ind,p.x,v);
  
  // outgoing velocity vector
  double beta = 0.0;
  for (int i=0;i<3;i++) 
  {
    v[i] = sign*v[i];
    beta += v[i]*v[i]; 
  }
  beta = sqrt(beta)/pc::c;
  double gamma  = 1.0/sqrt(1 - beta*beta); 

  // doppler shifts outgoing
  double vdp    = (p.D[0]*v[0] + p.D[1]*v[1] + p.D[2]*v[2]);
  double dshift = gamma*(1 - vdp/pc::c);

  // doppler shift the energy and frequency
  p.e   *= dshift;
  p.nu  *= dshift;

  // transform direction
  double D_old[3];
  D_old[0] = p.D[0];
  D_old[1] = p.D[1];
  D_old[2] = p.D[2];

  // See Mihalas & Mihalas eq 89.8
  p.D[0] = 1.0/dshift*(D_old[0]-gamma*v[0]/pc::c*(1-gamma*vdp/pc::c/(gamma+1)));
  p.D[1] = 1.0/dshift*(D_old[1]-gamma*v[1]/pc::c*(1-gamma*vdp/pc::c/(gamma+1)));
  p.D[2] = 1.0/dshift*(D_old[2]-gamma*v[2]/pc::c*(1-gamma*vdp/pc::c/(gamma+1)));

  // for security, make sure it is properly normalized
  double norm = p.D[0]*p.D[0] + p.D[1]*p.D[1] + p.D[2]*p.D[2];
  p.D[0] = p.D[0]/norm;
  p.D[1] = p.D[1]/norm;
  p.D[2] = p.D[2]/norm;
}


//----------------------------------------------------------------
//----------------------------------------------------------------
double species_general::sample_zone_nu(int zone_index)
{
  
}


//----------------------------------------------------------------
// return the emissivity integrated over nu for the core
//----------------------------------------------------------------
double species_general::int_core_emis()
{
  // TODO - implement photons::net_core_emis
}

//----------------------------------------------------------------
// return the emissivity integrated over nu for a zone
//----------------------------------------------------------------
double species_general::int_zone_emis(int zone_index)
{
  // TODO - implement photons::net_zone_emis
}


void species_general::propagate_particles(double dt)
{
  vector<particle>::iterator pIter = particles.begin();
  int n_active = particles.size();
  int n_escape = 0;
  while (pIter != particles.end())
    {
      ParticleFate fate = propagate(*pIter,dt);
      if (fate == escaped) n_escape++;
      if ((fate == escaped)||(fate == absorbed)) particles.erase(pIter);
      else pIter++;
    }
  double per_esc = (100.0*n_escape)/n_active;
  if (sim->verbose && sim->iterate) cout << "# Percent escaped = " << per_esc << endl;
  spectrum.rescale(1.0/per_esc);
}

//--------------------------------------------------------
// Propagate a single monte carlo particle until
// it  escapes, is absorbed, or the time step ends
//--------------------------------------------------------
ParticleFate species_general::propagate(particle &p, double dt)
{
  enum ParticleEvent {scatter, boundary, tstep};
  ParticleEvent event;

  ParticleFate  fate = moving;

  // time of end of timestep
  double tstop = t_now + dt;

  // local variables
  double tau_r,d_sc,d_tm,this_d;

  // pointer to current zone
  zone *zone = &(sim->grid->z[p.ind]);

  // propagate until this flag is set
  while (fate == moving)
  {
    // set pointer to current zone
    zone = &(sim->grid->z[p.ind]);
    
    // maximum step size inside zone
    double d_bn = sim->step_size * sim->grid->zone_min_length(p.ind);

    // doppler shift from comoving to lab
    double dshift = dshift_comoving_to_lab(p);

    // get local opacity and absorption fraction (epsilon)
    double e,a,s;
    get_eas(p,dshift,e,a,s);
    
    // convert opacity from comoving to lab frame for the purposes of 
    // determining the interaction distance in the lab frame
    // This corresponds to equation 90.8 in Mihalas&Mihalas. You multiply 
    // the comoving opacity by nu_0 over nu, which is why you
    // multiply by dshift instead of dividing by dshift here
    a = a*dshift;

    // random optical depth to next interaction
    tau_r = -1.0*log(1 - gsl_rng_uniform(sim->rangen));
    
    // step size to next interaction event
    d_sc  = tau_r/a;
    if (a == 0) d_sc = INFINITY;
    if (d_sc < 0) cout << "ERROR: negative interaction distance!\n" << endl;

    // find distance to end of time step
    d_tm = (tstop - p.t)*pc::c;
    if (sim->iterate) d_tm = INFINITY; // i.e. let all particles escape

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
    double this_E = p.e*this_d; 
    zone->e_rad += this_E; 

    // shift opacity back to comoving frame for energy and momentum exchange. 
    // Radiation energy is still lab frame
    a = a / dshift;

    // store absorbed energy in *comoving* frame 
    // (will turn into rate by dividing by dt later)
    // Extra dshift definitely needed here (two total)
    zone->e_abs  += this_E*dshift*(a)*eps*zone->eps_imc * dshift; 

    // put back in radiation force tally here
    // fx_rad =
    
    // move particle the distance
    p.x[0] += this_d*p.D[0];
    p.x[1] += this_d*p.D[1];
    p.x[2] += this_d*p.D[2]; 
    // advance the time
    p.t = p.t + this_d/pc::c;

    // ---------------------------------
    // Do if scatter
    // ---------------------------------
    if (event == scatter)
    {
      // random number to check for scattering or absorption
      double z = gsl_rng_uniform(sim->rangen);
      
      // do photon interaction physics
      if (p.type == photon)
      {
	// see if scattered 
	if (z > eps) isotropic_scatter(p,0);
	else
	{
	  // check for effective scattering
	  double z2;
	  if (sim->radiative_eq) z2 = 2;
	  else z2 = gsl_rng_uniform(sim->rangen);
	  // do an effective scatter
	  if (z2 > zone->eps_imc) isotropic_scatter(p,1);
	  // otherwise really absorb (kill) it
	  else fate = absorbed; 
	}
      }
    }
    
    // ---------------------------------
    // do if time step end
    // ---------------------------------
    else if (event == tstep) fate = stopped;

    // Find position of the particle now
    p.ind = sim->grid->get_zone(p.x);
    if (p.ind == -1) fate = absorbed;
    if (p.ind == -2) fate = escaped;

    // check for inner boundary absorption
    if (p.r() < sim->r_core) fate = absorbed;
  }

  // Add escaped photons to output spectrum
  if (fate == escaped) spectrum.count(p.t, p.nu, p.e, p.D); 

  return fate;
}

