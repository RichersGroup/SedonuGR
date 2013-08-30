#ifndef _SPECIES_H
#define _SPECIES_H

#include "particle.h"
#include "spectrum_array.h"
#include "Lua.h"
#include "grid_general.h"
#include "locate_array.h"
#include "cdf_array.h"
#include <gsl/gsl_rng.h>
#include <list>

class transport;

class species_general
{

 protected:

  // the frequency grid for emissivity/opacity
  locate_array nu_grid;

  // the core emissivity
  cdf_array core_emis;

  // the zone eas variables
  vector< cdf_array      > emis;
  vector< vector<double> > abs_opac;
  vector< vector<double> > scat_opac;

  // grey opacity and absorption fraction
  double grey_opac;
  double eps;

  // pointer to the simulation info (one level up)
  transport* sim;

 public:

  virtual void ~species_general();

  // name
  string name;

  // this species' spectrum
  spectrum_array spectrum;

  // lepton number (the particle property. really just 'electron number' used for
  // calculating changes to Ye. 1, -1, or 0)
  int lepton_number;

  // set everything up
  void init(Lua* lua, transport* sim);
  virtual void myInit(Lua* lua) = 0;

  // return the emissivity integrated over frequency at the core
  double int_core_emis();

  // return the emissivity integrated over frequency at a zone
  double int_zone_emis(int zone_index);
  double int_zone_lepton_emis(int zone_index);

  // return the frequency of a particle emitted from the core
  double sample_core_nu();

  // return the frequency of a particle emitted from a zone
  double sample_zone_nu(int zone_index);

  // set the emissivity, absorption opacity, and scattering opacity
  virtual void set_eas(int zone_index) = 0;
  void get_opacity(particle &p, double dshift, double* opac, double* abs_frac);

  // min and max values for the Brent solver
  double T_min,  T_max;
  double Ye_min, Ye_max;
};




#endif
