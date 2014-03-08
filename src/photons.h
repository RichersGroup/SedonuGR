#ifndef _PHOTONS_H
#define _PHOTONS_H

#include "species_general.h"
#include "Lua.h"

class photons: public species_general
{

protected:

  // photon-specific functions
  double klein_nishina(const double) const;
  double planck(const double T, const double nu) const; // (erg/s/cm^2/Hz/ster)
  void compton_scatter();

public:

  virtual ~photons() {}

  // required functions
  void myInit(Lua* lua);
  void set_eas(int zone_index);
};

#endif
