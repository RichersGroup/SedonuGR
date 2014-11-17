#ifndef _SPECIES_H
#define _SPECIES_H

#include "global_options.h"
#include <list>
#include "particle.h"
#include "spectrum_array.h"
#include "Lua.h"
#include "grid_general.h"
#include "locate_array.h"
#include "cdf_array.h"

class transport;

class species_general
{

protected:

	// the frequency grid for emissivity/opacity (Hz)
	locate_array nu_grid;

	// the zone eas variables
	vector< cdf_array      > emis;
	vector< vector<double> > abs_opac;  // 1/cm
	vector< vector<double> > scat_opac; // 1/cm

	// grey opacity and absorption fraction
	double grey_opac; //(cm^2/g)
	double grey_abs_frac;       //unitless

	// pointer to the simulation info (one level up)
	transport* sim;

	// species-specific initialization stuff
	virtual void myInit(Lua* lua) = 0;

public:

	species_general();
	virtual ~species_general() {}

	// name
	string name;

	// lepton number (the particle property, {-1,0,1})
	int lepton_number;

	// the numbers of species this represents
	double weight;

	// this species' spectrum
	spectrum_array spectrum;

	// the core emissivity (units of core_emis.N are erg/s)
	cdf_array core_emis;

	// set everything up
	void init(Lua* lua, transport* sim);

	// this species' blackbody function (erg/cm^2/s/ster/Hz)
	virtual double blackbody(const double T, const double chempot, const double nu) const = 0;

	// set a CDF to blackbody distribution
	void set_cdf_to_BB(const double T, const double chempot, cdf_array& emis);

	// return the emissivity integrated over frequency at the core
	double integrate_core_emis() const; //(erg/s)

	// return the emissivity integrated over frequency at a zone
	double integrate_zone_emis(const int zone_index) const;        //(erg/s/cm^3/ster)
	double integrate_zone_lepton_emis(const int zone_index) const; //unitless

	// return the frequency of a particle emitted from the core (Hz)
	double sample_core_nu(const int g=-1) const;

	// return the frequency of a particle emitted from a zone (Hz)
	double sample_zone_nu(const int zone_index, const int g=-1) const;

	// set the emissivity, absorption opacity, and scattering opacity
	virtual void set_eas(const int zone_index) = 0;
	void set_emis_to_BB_edens(const double T, const double chempot);
	void get_opacity(const double com_nu, const int z_ind, double* opac, double* abs_frac) const;
	double sum_opacity(const int z_ind, const int group) const;

	// minimum zone emissivity
	double bin_emis(const int zone_index, const int g) const;
	double min_bin_emis(const int zone_index) const;
	unsigned number_of_bins();

	// min and max values for the Brent solver
	double T_min,  T_max; //(K)
	double Ye_min, Ye_max;
	double rho_min, rho_max; //(g/cm^3)
};




#endif
