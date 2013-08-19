#ifndef _NEUTRINOS_H
#define _NEUTRINOS_H

#include "species_general.h"
#include "locate_array.h"
#include "cdf_array.h"
#include <vector>
#include "Lua.h"
#include <gsl/gsl_rng.h>

class neutrinos: public species_general
{

 protected:  

 public:

  int num_nut_species;
  int nulibID;
  int electron_number;

  // required functions
  void myInit(Lua* lua);
  void set_eas(int zone_index);
};

#endif