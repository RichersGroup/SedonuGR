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

#ifndef _NEUTRINOS_H
#define _NEUTRINOS_H

#include "Species.h"

// PARAMETERS
//    Neutrino_cdf_cutoff - cutoff bin value in CDF. Below this value is considered zero. Necessary for some equilibrium problems.
//    Neutrino_grey_opac - if <0, takes opac/emis from NuLib. Else uses this constant opacity.
//    Neutrino_grey_abs_frac - (if grey_opac<0) constant absorption fraction
//    Neutrino_nugrid_start - (if grey_opac<0) bottom of frequency grid
//    Neutrino_nugrid_stop - (if grey_opac<0) top of freqyency grid
//    Neutrino_nugrid_n - (if grey_opac<0) number of frequency bins
//    Neutrino_spec_n_mu - # of mu bins in neutrino escape spectrum
//    Neutrino_spec_n_phi - # of phi bins in neutrino escape spectrum

class Neutrino: public Species
{

protected:

public:

	Neutrino();
	virtual ~Neutrino() {}

	int num_species;

	// virtual functions
	virtual void set_eas(int zone_index) = 0;

	// required functions
	void myInit(Lua* lua);
	double blackbody(const double T, const double chempot, const double nu) const;
	static double annihilation_rate(const SpectrumArray& nu_dist, const SpectrumArray& nbar_dist, const bool electron_type, const int weight);
};

#endif
