/*
//  Copyright (c) 2015, California Institute of Technology and the Regents
//  of the University of California, based on research sponsored by the
//  United States Department of Energy. All rights reserved.
//
//  This file is part of Sedonu.
//
//  Sedonu is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  Neither the name of the California Institute of Technology (Caltech)
//  nor the University of California nor the names of its contributors 
//  may be used to endorse or promote products derived from this software
//  without specific prior written permission.
//
//  Sedonu is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with Sedonu.  If not, see <http://www.gnu.org/licenses/>.
//
*/

#ifndef _SPECIES_H
#define _SPECIES_H

#include <vector>
#include "LorentzHelper.h"
#include "SpectrumArray.h"
#include "Lua.h"
#include "CDFArray.h"

class Transport;

class Species
{

protected:

	// the frequency grid for emissivity/opacity (Hz)
	LocateArray nu_grid;

	// the zone eas variables
	std::vector< CDFArray      > emis;
	std::vector< CDFArray      > biased_emis;
	std::vector< std::vector<double> > abs_opac;  // 1/cm
	std::vector< std::vector<double> > scat_opac; // 1/cm

	// pointer to the simulation info (one level up)
	Transport* sim;

	// species-specific initialization stuff
	virtual void myInit(Lua* lua) = 0;

public:

	Species();
	virtual ~Species() {}

	// name
	std::string name;

	// lepton number (the particle property, {-1,0,1})
	int lepton_number;

	// the numbers of species this represents
	double weight;

	// this species' spectrum
	SpectrumArray spectrum;

	// the core emissivity (units of core_emis.N are erg/s)
	CDFArray core_emis;

	// set everything up
	void init(Lua* lua, Transport* sim);

	// this species' blackbody function (erg/cm^2/s/ster/Hz)
	virtual double blackbody(const double T, const double chempot, const double nu) const = 0;

	// set a CDF to blackbody distribution
	void set_cdf_to_BB(const double T, const double chempot, CDFArray& emis);

	// return the emissivity integrated over frequency at the core
	double integrate_core_emis() const; //(erg/s)

	// return the emissivity integrated over frequency at a zone
	double integrate_zone_emis(const int zone_index) const;        //(erg/s/cm^3/ster)
	double integrate_zone_lepton_emis(const int zone_index) const; //unitless
	double integrate_zone_biased_emis(const int zone_index) const; //(erg/s/cm^3/ster)

	// return the frequency of a particle emitted from the core or a zone (Hz)
	double sample_core_nu(const int g=-1) const;
	double sample_nu(const CDFArray& input_emis, const int g=-1) const;
	double sample_zone_nu(const int zone_index, double *Ep, const int g=-1) const;

	// set the emissivity, absorption opacity, and scattering opacity
	virtual void set_eas(const int zone_index) = 0;
	void get_opacity(const double com_nu, const int z_ind, double* abs_opac, double* scat_opac) const;
	double sum_opacity(const int z_ind, const int group) const;
	double interpolate_importance(const double nu, const int z_ind) const;

	// minimum zone emissivity
	double bin_emis(const int zone_index, const int g) const;
	unsigned number_of_bins();

	// min and max values for the Brent solver
	double T_min,  T_max; //(K)
	double Ye_min, Ye_max;
	double rho_min, rho_max; //(g/cm^3)
};




#endif
