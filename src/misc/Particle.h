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

#ifndef _PARTICLE_H
#define _PARTICLE_H

#include <cmath>
#include <cstdio>
#include "global_options.h"

enum ParticleFate  {moving, escaped, absorbed, rouletted};

// particle class
class Particle
{

public:

	Particle(){
		e = NaN;
		nu = NaN;
		tau = NaN;
		s = -1;
		fate = moving;
	}

	double x[3];            // x,y,z position (cm)
	double D[3];            // direction vector, Dx,Dy,Dz (normalized to unit magnitude)
	double       e;         // total energy in ergs of packet
	double      nu;         // frequency (Hz)
	double     tau;         // remaining optical depth
	int          s;         // species number
	ParticleFate fate;

	void print() const
	{
		printf("%10.3e %10.3e %10.3e | %10.3e %10.3e %10.3e | %10.3e %10.3e %10.3e %i %i\n",
				x[0],x[1],x[2],D[0],D[1],D[2],e,nu,tau,s,fate);
	}

};

#endif
