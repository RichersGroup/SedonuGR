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
#include <iostream>
#include "global_options.h"

enum ParticleFate  {moving, escaped, absorbed, rouletted};

// particle class
class Particle
{

public:

	Particle(){
		for(int i=0; i<4; i++){
			xup[i] = NaN;
			kup[i] = NaN;
		}
		N = NaN;
		tau = NaN;
		s = -1;
		fate = moving;
	}

	double xup[4];          // x,y,z,ct position (cm)
	double kup[4];          // 4-wavevector (2pi nu/c)
	double       N;         // total number of neutrinos in packet
	double     tau;         // remaining optical depth
	int          s;         // species number
	ParticleFate fate;

	void print() const
	{
		printf("%10.3e %10.3e %10.3e %10.3e | %10.3e %10.3e %10.3e %10.3e | %10.3e %10.3e %i %i\n",
				xup[0],xup[1],xup[2],xup[3],kup[0],kup[1],kup[2],kup[3],N,tau,s,fate);
	}

};

#endif
