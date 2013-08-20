#include <math.h>
#include <gsl/gsl_rng.h>
#include "transport.h"
#include "species_general.h"
#include "physical_constants.h"

namespace pc = physical_constants;

//------------------------------------------------------------
// emit new particles
//------------------------------------------------------------
void transport::emit_particles(double dt)
{
  if(n_inject > 0) emit_inner_source(dt);
  // TODO - implement emitting from zones when not doing radiative equilibrium
}


//------------------------------------------------------------
// General function to create a particle in zone i
// emitted isotropically in the comoving frame. 
// Useful for thermal radiation emitted all througout
// the grid
//------------------------------------------------------------
void transport::create_isotropic_particle(int zone_index, double Ep)
{
  // basic particle properties
  particle p;
  p.ind = zone_index;
  p.e  = Ep;
  p.t  = t_now;

  // random sample position in zone
  vector<double> rand(3,0);
  rand[0] = gsl_rng_uniform(rangen);
  rand[1] = gsl_rng_uniform(rangen);
  rand[2] = gsl_rng_uniform(rangen);
  double r[3];
  grid->sample_in_zone(zone_index,rand,r);
  p.x[0] = r[0];
  p.x[1] = r[1];
  p.x[2] = r[2];

  // emit isotropically in comoving frame
  double mu  = 1 - 2.0*gsl_rng_uniform(rangen);
  double phi = 2.0*pc::pi*gsl_rng_uniform(rangen);
  double smu = sqrt(1 - mu*mu);
  p.D[0] = smu*cos(phi);
  p.D[1] = smu*sin(phi);
  p.D[2] = mu;

  // sample the species and frequency
  int s = sample_zone_species(zone_index);
  p.nu = species_list[s]->sample_zone_nu(zone_index);

  // subtract the leptons from the zone
  grid->z[p.ind].l_abs -= species_list[s]->lepton_number * p.e/(p.nu*pc::h);

  // lorentz transform from the comoving to lab frame
  species_list[s]->transform_comoving_to_lab(p);

  // add to particle vector
  species_list[s]->add_particle(p);
}


//------------------------------------------------------------
// Initialize a constant number of particles
// in each zone
//------------------------------------------------------------
void transport::initialize_particles(int init_particles)
{
  if (verbose) cout << "# initializing with " << init_particles << " particle per zone\n";

  for (int i=0;i<grid->z.size();i++)
  {
    // lab frame energy
    double E_zone = grid->z[i].e_rad*grid->zone_volume(i);
    // particle energy
    double Ep = E_zone/(init_particles);

    // create init_particles particles
    for (int q=0;q<init_particles;q++) create_isotropic_particle(i,Ep);
  }
}


//------------------------------------------------------------
// inject particles from a central luminous source
// Currently written to emit photons with 
// blackblody spectrum based on T_core and L_core
//------------------------------------------------------------
void transport::emit_inner_source(double dt)
{
  if (total_particles() + n_inject > MAX_PARTICLES){
    cout << "# Not enough particle space\n"; 
    exit(10);
  }

  double Ep  = L_core*dt/n_inject;
  int s;

  // inject particles from the source
  for (int i=0;i<n_inject;i++)
  {
    // set basic properties
    particle p;
    
    // pick initial position on photosphere
    double phi_core   = 2*pc::pi*gsl_rng_uniform(rangen);
    double cosp_core  = cos(phi_core);
    double sinp_core  = sin(phi_core);
    double cost_core  = 1 - 2.0*gsl_rng_uniform(rangen);
    double sint_core  = sqrt(1-cost_core*cost_core);
    // real spatial coordinates    
    double a_phot = r_core + r_core*1e-10;
    p.x[0] = a_phot*sint_core*cosp_core;
    p.x[1] = a_phot*sint_core*sinp_core;
    p.x[2] = a_phot*cost_core;

    // pick photon propagation direction wtr to local normal                   
    double phi_loc = 2*pc::pi*gsl_rng_uniform(rangen);
    // choose sqrt(R) to get outward, cos(theta) emission         
    double cost_loc  = sqrt(gsl_rng_uniform(rangen));
    double sint_loc  = sqrt(1 - cost_loc*cost_loc);
    // local direction vector                     
    double D_xl = sint_loc*cos(phi_loc);
    double D_yl = sint_loc*sin(phi_loc);
    double D_zl = cost_loc;
    // apply rotation matrix to convert D vector into overall frame        
    p.D[0] = cost_core*cosp_core*D_xl-sinp_core*D_yl+sint_core*cosp_core*D_zl;
    p.D[1] = cost_core*sinp_core*D_xl+cosp_core*D_yl+sint_core*sinp_core*D_zl;
    p.D[2] = -sint_core*D_xl+cost_core*D_zl;

    // set basic properties
    p.e = Ep;
    p.t  = t_now;

    // get index of current zone
    p.ind = grid->get_zone(p.x);
    if(p.ind < 0){
      cout << "WARNING: particle spawned with emit_inner_source is outside the grid" << endl;
    }

    // sample the species and frequency
    int s = sample_core_species();
    p.nu = species_list[s]->sample_core_nu();

    // lorentz transform from the comoving to lab frame
    species_list[s]->transform_comoving_to_lab(p);

    // add to particle vector
    species_list[s]->add_particle(p);
  }

  if (verbose) printf("# Injected %d particles\n",n_inject);
}


