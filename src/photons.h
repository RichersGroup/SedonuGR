#ifndef _PHOTONS_H
#define _PHOTONS_H

#include "species_general.h"
#include "Lua.h"

class photons: public species_general
{

protected:

  // photon-specific functions
  double klein_nishina(const double) const;
  void compton_scatter();

public:

  virtual ~photons() {}

  // required functions
  void myInit(Lua* lua);
  void set_eas(int zone_index);
  double blackbody(const double T, const double chempot, const double nu) const;
};

#endif
