#pragma warning disable 161
#include <omp.h>
#include <math.h>
#include "thread_RNG.h"
#include "transport.h"
#include "species_general.h"
#include "physical_constants.h"

namespace pc = physical_constants;

//------------------------------------------------------------
// Initialize a constant number of particles in each zone
//------------------------------------------------------------
void transport::initialize_particles(const int init_particles){
  if (verbose) cout << "# initializing with " << init_particles << " particle per zone\n";

  #pragma omp parallel for
  for (int i=0;i<grid->z.size();i++){
    double E_zone = grid->z[i].e_rad*grid->zone_volume(i);  // lab frame energy
    double Ep = E_zone/init_particles;                      // particle energy
    for (int q=0;q<init_particles;q++) create_thermal_particle(i,Ep,t_now);
  }
}


//------------------------------------------------------------
// emit new particles
//------------------------------------------------------------
void transport::emit_particles(const double dt)
{
  int old_size = particles.size();
  L_net = 0; // re-summed in emission routines

  // complain if we're out of room for particles
  int n_emit = n_emit_core + n_emit_therm + n_emit_decay;
  if (total_particles() + n_emit > max_particles){
    cout << "# ERROR: Not enough particle space\n";
    exit(10);
  }

  // emit from the core and/or the zones
  if(do_core) emit_inner_source(dt);
  if(do_therm || do_decay) emit_zones(dt);

  // print how many particles were added to the system
  int new_size = particles.size();
  if (verbose) printf("# Emitted %d particles\n",new_size-old_size);
}


//------------------------------------------------------------
// inject particles from a central luminous source
// Currently written to emit photons with 
// blackblody spectrum based on T_core and L_core
//------------------------------------------------------------
void transport::emit_inner_source(const double dt)
{
  const double Ep  = L_core*dt/n_emit_core;
  L_net += L_core;

  #pragma omp parallel for
  for (int i=0; i<n_emit_core; i++){
    double t = t_now + dt*rangen.uniform();
    create_surface_particle(Ep,t);
  }
}


//--------------------------------------------------------------------------
// emit particles from the zones due to non-radiation heating
//--------------------------------------------------------------------------
void transport::emit_zones(const double dt)
{
  int gridsize = grid->z.size(); 
  double therm_lum = 0;
  double decay_lum = 0;
  double Ep_decay, Ep_therm;
  int asdfindex = 230400;

  // at this point therm means either viscous heating or regular emission, according to the logic above
  #pragma omp parallel
  {

    // determine the net luminosity of each emission type over the whole grid
    #pragma omp for reduction(+:therm_lum,decay_lum)
    for(int i=0; i<gridsize; i++){
      if(i==asdfindex){
	cout << (radiative_eq ? zone_visc_heat_rate(asdfindex) : zone_heat_lum(asdfindex)) << endl;// TODO delete test
	cout << (radiative_eq ? zone_visc_heat_rate(asdfindex) : 0 /*zone_heat_lum(asdfindex)*/) << endl;// TODO delete test
	cout << (radiative_eq ? zone_visc_heat_rate(asdfindex) : zone_heat_lum(asdfindex)) << endl;// TODO delete test
      }
      if(do_therm) therm_lum += ( radiative_eq ? zone_visc_heat_rate(i) : zone_heat_lum(i) );
      if(do_decay) decay_lum += zone_decay_lum(i);
    }

    // store the variables in transport and set the particle energy
    #pragma omp single
    {
      if(do_therm && !radiative_eq) cout << "L_therm = " << therm_lum << endl;
      L_net += therm_lum + decay_lum;
      if(do_therm) Ep_therm = therm_lum*dt / (double)n_emit_therm;
      if(do_decay) Ep_decay = decay_lum*dt / (double)n_emit_decay;
      cout << (radiative_eq ? zone_visc_heat_rate(asdfindex) : zone_heat_lum(asdfindex)) << endl;// TODO delete test
    }
    
    // create particles in each grid cell
    #pragma omp for schedule(guided)
    double asdf=0;
    for (int i=0; i<gridsize; i++)
    {

      // EMIT THERMAL PARTICLES =========================================================
      if(do_therm && therm_lum>0){
	// this zone's luminosity and number of emitted particles.
	// randomly decide whether last particle gets added based on the remainder.
	if(i==asdfindex){
	  cout << (radiative_eq ? zone_visc_heat_rate(i) : zone_heat_lum(i)) << endl;// TODO delete test
	}
	double this_L  = ( radiative_eq ? zone_visc_heat_rate(i) : zone_heat_lum(i) );
	double tmp  = (double)n_emit_therm * this_L/therm_lum;
	asdf += fmod(tmp,1.0);
	int this_n_emit = (int)tmp + (int)( rangen.uniform() < fmod(tmp,1.0) );

	// add heat absorbed to tally of e_abs
	// really, this is "the gas absorbs energy from heating, then emits radiation"
	// hence, this is only done if we assume radiative equilibrium
	if(radiative_eq) grid->z[i].e_abs += dt * this_L;

	// create the particles
	double t;
	for (int k=0; k<this_n_emit; k++){
	  t = t_now + dt*rangen.uniform();
	  create_thermal_particle(i,Ep_therm,t);
	}
      }

      // EMIT DECAY PARTICLES ===========================================================
      if(do_decay && decay_lum>0){
	// this zone's luminosity and number of emitted particles
	double this_L = zone_decay_lum(i);
	double tmp  = (double)n_emit_decay * this_L/decay_lum;
	int this_n_emit = (int)tmp + (int)( rangen.uniform() < fmod(tmp,1.0) );

	// create the particles
	double t;
	for (int k=0; k<this_n_emit; k++){
	  t = t_now + dt*rangen.uniform();
	  create_decay_particle(i,Ep_decay,t);
	}
      }

    }//loop over zones  
    cout << "asdf = " << asdf << endl;
  }//#pragma omp parallel
}


//----------------------------------------------------------------------------------------
// Helper functions for emit_zones
//----------------------------------------------------------------------------------------

// rate at which viscosity energizes the fluid (erg/s)
double transport::zone_visc_heat_rate(const int zone_index) const{
  return visc_specific_heat_rate * grid->z[zone_index].rho * grid->zone_volume(zone_index);
}

// return the cell's luminosity from thermal emission (erg/s)
double transport::zone_heat_lum(const int zone_index) const{
  double H;
  for(int i=0; i<species_list.size(); i++)
    H += species_list[i]->int_zone_emis(zone_index) * 4*pc::pi * grid->zone_volume(zone_index);
  return H;
}

// return the total decay luminosity from the cell (erg/s)
double transport::zone_decay_lum(const int zone_index) const{
  if(verbose) cout << "WARNING: emission in zones from decay is not yet implemented!" << endl;
  return 0;
}


//------------------------------------------------------------
// General function to create a particle in zone i
// emitted isotropically in the comoving frame. 
// Useful for thermal radiation emitted all througout
// the grid
//------------------------------------------------------------
void transport::create_thermal_particle(const int zone_index, const double Ep, const double t)
{
  // basic particle properties
  particle p;
  p.ind = zone_index;
  p.e  = Ep;
  p.t  = t;

  // random sample position in zone
  vector<double> rand(3,0);
  rand[0] = rangen.uniform();
  rand[1] = rangen.uniform();
  rand[2] = rangen.uniform();
  double r[3];
  grid->sample_in_zone(zone_index,rand,r);
  p.x[0] = r[0];
  p.x[1] = r[1];
  p.x[2] = r[2];

  // emit isotropically in comoving frame
  double mu  = 1 - 2.0*rangen.uniform();
  double phi = 2.0*pc::pi*rangen.uniform();
  double smu = sqrt(1 - mu*mu);
  p.D[0] = smu*cos(phi);
  p.D[1] = smu*sin(phi);
  p.D[2] = mu;

  // sample the species and frequency
  int s = sample_zone_species(zone_index);
  p.s = s;
  p.nu = species_list[s]->sample_zone_nu(zone_index);

  // lorentz transform from the comoving to lab frame
  transform_comoving_to_lab(&p);

  // add to particle vector
  #pragma omp critical
  particles.push_back(p);
}


//------------------------------------------------------------
// General function to create a particle on the surface
// emitted isotropically outward in the comoving frame. 
//------------------------------------------------------------
void transport::create_surface_particle(const double Ep, const double t)
{
  // set basic properties
  particle p;
  p.e = Ep;
  p.t = t;

  // pick initial position on photosphere
  double phi_core   = 2*pc::pi*rangen.uniform();
  double cosp_core  = cos(phi_core);
  double sinp_core  = sin(phi_core);
  double cost_core  = 1 - 2.0*rangen.uniform();
  double sint_core  = sqrt(1-cost_core*cost_core);
  // real spatial coordinates    
  double a_phot = r_core + r_core*1e-10;
  p.x[0] = a_phot*sint_core*cosp_core;
  p.x[1] = a_phot*sint_core*sinp_core;
  p.x[2] = a_phot*cost_core;

  // pick photon propagation direction wtr to local normal                   
  double phi_loc = 2*pc::pi*rangen.uniform();
  // choose sqrt(R) to get outward, cos(theta) emission         
  double cost_loc  = sqrt(rangen.uniform());
  double sint_loc  = sqrt(1 - cost_loc*cost_loc);
  // local direction vector                     
  double D_xl = sint_loc*cos(phi_loc);
  double D_yl = sint_loc*sin(phi_loc);
  double D_zl = cost_loc;
  // apply rotation matrix to convert D vector into overall frame        
  p.D[0] = cost_core*cosp_core*D_xl-sinp_core*D_yl+sint_core*cosp_core*D_zl;
  p.D[1] = cost_core*sinp_core*D_xl+cosp_core*D_yl+sint_core*sinp_core*D_zl;
  p.D[2] = -sint_core*D_xl+cost_core*D_zl;

  // get index of current zone
  p.ind = grid->get_zone(p.x);
  if(p.ind < 0){
    printf("WARNING: particle spawned with emit_inner_source is outside the grid.\n");
  }

  // sample the species and frequency
  int s = sample_core_species();
  p.s = s;
  p.nu = species_list[s]->sample_core_nu();

  // lorentz transform from the comoving to lab frame
  transform_comoving_to_lab(&p);

  // add to particle vector
  #pragma omp critical
  particles.push_back(p);
}


//------------------------------------------------------------
// General function to create a particle from radioactive decay
//------------------------------------------------------------
void transport::create_decay_particle(const int zone_index, const double Ep, const double t)
{
  if(verbose) cout << "WARNING: transport::create_decay_particle is not yet implemented!" << endl;
}
