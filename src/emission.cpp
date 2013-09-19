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
void transport::initialize_particles(int init_particles){
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
void transport::emit_particles(double dt)
{
  // complain if we're out of room for particles
  int n_emit = n_emit_core + n_emit_heat + n_emit_decay;
  if (total_particles() + n_emit > MAX_PARTICLES){
    cout << "# Not enough particle space\n"; 
    exit(10);
  }

  // set the net luminosities of the thermal and decay particles
  // if not iterative, the problem is time-varying, so the luminosity must be recalculated.
  // if not radiative_eq, the temperature determines the emissivity, so the luminosity must be recalculated.
  if(!iterate && !radiative_eq){
    double heat_lum = 0;
    double decay_lum = 0;
    #pragma omp parallel for reduction(+:heat_lum,decay_lum)
    for(int i=0; i<grid->z.size(); i++){
      L_heat += zone_heating_rate(i);
      L_decay += zone_decay_rate(i);
    }
    L_heat = heat_lum;
    L_decay = decay_lum;
  }

  // emit from the core and/or the zones
  if(n_emit_core > 0) emit_inner_source(dt);
  if(n_emit_heat>0 || n_emit_decay>0) emit_zones(dt);

  // print how many particles were added to the system
  if (verbose) printf("# Emitted %d particles\n",n_emit);
}


//------------------------------------------------------------
// inject particles from a central luminous source
// Currently written to emit photons with 
// blackblody spectrum based on T_core and L_core
//------------------------------------------------------------
void transport::emit_inner_source(double dt)
{
  const double Ep  = L_core*dt/n_emit_core;

  #pragma omp parallel for
  for (int i=0; i<n_emit_core; i++){
    double t = t_now + dt*rangen.uniform();
    create_surface_particle(Ep,t);
  }
}


//------------------------------------------------------------
// emit particles from the zones due to non-radiation heating
//------------------------------------------------------------
void transport::emit_zones(double dt)
{
  double t;
  double Ep_heat  = L_heat*dt / n_emit_heat;
  double Ep_decay = L_decay*dt / n_emit_decay;
  
  #pragma omp parallel for schedule(guided)
  for (int i=0;i<grid->z.size();i++)
  {
    // EMIT THERMAL PARTICLES =========================================================
    if(L_heat>0 && n_emit_heat>0){
      // this zone's luminosity and number of emitted particles
      double this_L_heat = zone_heating_rate(i);
      int this_n_emit = n_emit_heat * (this_L_heat/L_heat + 0.5);

      // add to tally of e_abs
      // really, this is "the gas absorbs energy from heating, then emits radiation"
      // hence, this is only done if we assume radiative equilibrium
      if(radiative_eq) grid->z[i].e_abs += dt * this_L_heat;
      
      // create the particles
      for (int k=0; k<this_n_emit; k++){
	double t = t_now + dt*rangen.uniform();
	create_thermal_particle(i,Ep_heat,t);
      }
    }


    // EMIT DECAY PARTICLES ===========================================================
    if(L_decay>0 && n_emit_decay>0){
      // this zone's luminosity and number of emitted particles
      double this_L_decay = zone_decay_rate(i);
      int this_n_emit = n_emit_decay * (this_L_decay/L_decay + 0.5);

      // create the particles
      for (int k=0; k<this_n_emit; k++){
	double t = t_now + dt*rangen.uniform();
	create_decay_particle(i,Ep_decay,t);
      }
    }

  } // loop over zones  
}


// return the total heating rate in the cell in erg/s
double transport::zone_heating_rate(int zone_index)
{
  if(radiative_eq) return 0;//specific_heating_rate * grid->z[zone_index].rho * grid->zone_volume(zone_index);
  else{
    double H;
    for(int i=0; i<species_list.size(); i++)
      H += species_list[i]->int_zone_emis(zone_index) * 4*pc::pi * grid->zone_volume(zone_index);
    return H;
  }
}


// return the total decay luminosity from the cell in erg/s
double transport::zone_decay_rate(int zone_index)
{
  if(verbose) cout << "WARNING: emission in zones from decay is not yet implemented!" << endl;
  return 0;
}


//------------------------------------------------------------
// General function to create a particle in zone i
// emitted isotropically in the comoving frame. 
// Useful for thermal radiation emitted all througout
// the grid
//------------------------------------------------------------
void transport::create_thermal_particle(int zone_index, double Ep, double t)
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
void transport::create_surface_particle(double Ep, double t)
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
void transport::create_decay_particle(int zone_index, double Ep, double t)
{
  if(verbose) cout << "WARNING: transport::create_decay_particle is not yet implemented!" << endl;
}
