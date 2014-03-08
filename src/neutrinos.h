#ifndef _NEUTRINOS_H
#define _NEUTRINOS_H

#include "species_general.h"
#include "Lua.h"

class neutrinos: public species_general
{

 protected:  

 public:

  virtual ~neutrinos() {}

  int num_nut_species;
  int nulibID;

  // required functions
  void myInit(Lua* lua);
  void set_eas(int zone_index);

  // other functions
  double fermi_dirac(const double T, const double chem_pot, const double nu) const; //(unitless)
};

#endif
