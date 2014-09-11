#ifndef _PARTICLE_H
#define _PARTICLE_H
#include <math.h>
#include <stdio.h>
#include <vector>
#include "global_options.h"

enum ParticleFate  {moving, stopped, escaped, absorbed};

// particle class
class particle
{

public:

	particle();

	vector<double> x;         // x,y,z position (cm)
	vector<double> D;         // direction vector, Dx,Dy,Dz (normalized to unit magnitude)
	double       t;         // current time
	double       e;         // total energy in ergs of packet
	double      nu;         // frequency (Hz)
	int          s;         // species number
	ParticleFate fate;

	double r() const
	{ return sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]); }

	double x_dot_d() const
	{return x[0]*D[0] + x[1]*D[1] + x[2]*D[2]; }

	double mu() const{
		double radius = r();
		if(radius == 0) return 0;
		else return x_dot_d()/r();
	}

	void normalize_direction();

	void print() const
	{
		printf("%10.e | %10.3e %10.3e %10.3e | %10.3e %10.3e %10.3e | %10.3e %10.3e\n",
				t,x[0],x[1],x[2],D[0],D[1],D[2],e,nu);
	}

};

#endif
