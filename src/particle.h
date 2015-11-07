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
#include <math.h>
#include <stdio.h>
#include <vector>
#include "global_options.h"

enum ParticleFate  {moving, escaped, absorbed};

// particle class
class particle
{

public:

	particle();

	vector<double> x;         // x,y,z position (cm)
	vector<double> D;         // direction vector, Dx,Dy,Dz (normalized to unit magnitude)
	double       e;         // total energy in ergs of packet
	double      nu;         // frequency (Hz)
	double     tau;         // remaining optical depth
	int          s;         // species number
	ParticleFate fate;

	double r() const
	{ return sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]); }

	double rcyl() const
	{ return sqrt(x[0]*x[0] + x[1]*x[1]); }

	double x_dot_d() const
	{return x[0]*D[0] + x[1]*D[1] + x[2]*D[2]; }

	double xcyl_dot_dcyl() const
	{return x[0]*D[0] + x[1]*D[1]; }

	double mu() const{
		double radius = r();
		if(radius == 0) return 0;
		else return x_dot_d()/r();
	}

	double mucyl() const{
		double radius = rcyl();
		if(radius == 0) return 0;
		else return xcyl_dot_dcyl()/rcyl();
	}

	void print() const
	{
		printf("%10.3e %10.3e %10.3e | %10.3e %10.3e %10.3e | %10.3e %10.3e\n",
				x[0],x[1],x[2],D[0],D[1],D[2],e,nu);
	}

};

#endif
