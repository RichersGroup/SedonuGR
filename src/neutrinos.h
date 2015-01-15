#ifndef _NEUTRINOS_H
#define _NEUTRINOS_H

#include "species_general.h"
#include "Lua.h"
#include "global_options.h"

class neutrinos: public species_general
{

protected:

public:

	neutrinos();
	virtual ~neutrinos() {}

	int num_nut_species;
	int nulibID;
	double cutoff;

	// required functions
	void myInit(Lua* lua);
	void set_eas(int zone_index);
	double blackbody(const double T, const double chempot, const double nu) const;
	static double annihilation_rate(const spectrum_array& nu_dist, const spectrum_array& nbar_dist, const bool electron_type, const int weight);
};

#endif
