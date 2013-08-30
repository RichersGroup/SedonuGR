#ifndef _PHOTONS_H
#define _PHOTONS_H

#include "species_general.h"
#include "locate_array.h"
#include "cdf_array.h"
#include <vector>
#include "Lua.h"
#include <gsl/gsl_rng.h>

class photons: public species_general
{

protected:

  // photon-specific functions
  double klein_nishina(double);
  double planck(double T, double nu);
  void compton_scatter();

public:

  virtual ~photons() {}

  // required functions
  void myInit(Lua* lua);
  void set_eas(int zone_index);
};

#endif
