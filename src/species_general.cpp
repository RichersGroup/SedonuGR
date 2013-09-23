#include "species_general.h"
#include <gsl/gsl_rng.h>
#include <math.h>
#include "transport.h"
#include "physical_constants.h"
#include "Lua.h"
#include <list>

namespace pc = physical_constants;

void species_general::init(Lua* lua, transport* simulation)
{
  // set the pointer to see the simulation info
  sim = simulation;

  // call child's init function
  myInit(lua);
}




//----------------------------------------------------------------
// return a randomly sampled frequency
// for a particle emitted from the core
//----------------------------------------------------------------
// TODO - could be modified to interpolate via inverse transform sampling, but more complicated
double species_general::sample_core_nu()
{
  // randomly pick a frequency bin
  double z = sim->rangen.uniform();
  int ilam = core_emis.sample(z);

  // randomly pick a location in the frequency bin
  z = sim->rangen.uniform();
  return nu_grid.sample(ilam,z);
}

//----------------------------------------------------------------
// return a randomly sampled frequency
// for a particle emitted from a zone
//----------------------------------------------------------------
// TODO - could be modified to interpolate via inverse transform sampling, but more complicated
double species_general::sample_zone_nu(int zone_index)
{
  // randomly pick a frequency bin
  double z = sim->rangen.uniform();
  int ilam = emis[zone_index].sample(z);

  // randomly pick a location in the frequency bin
  z = sim->rangen.uniform();
  return nu_grid.sample(ilam,z);
}


//----------------------------------------------------------------
// return the emissivity integrated over nu for the core (erg/s)
//----------------------------------------------------------------
double species_general::int_core_emis()
{
  return core_emis.N;
}

//----------------------------------------------------------------
// return the emissivity integrated over nu for a zone (erg/s/ster/cm^3)
//----------------------------------------------------------------
double species_general::int_zone_emis(int zone_index)
{
  return emis[zone_index].N;
}


//----------------------------------------------------------------
// return the lepton emissivity integrated over nu for a zone (#/s/ster/cm^3)
//----------------------------------------------------------------
double species_general::int_zone_lepton_emis(int zone_index)
{
  double l_emis = 0;
  for(int i=0; i<emis[zone_index].size(); i++)
  {
    l_emis += lepton_number * emis[zone_index].get_value(i) / (pc::h*nu_grid.x[i]);
  }
  return l_emis * emis[zone_index].N;
}
