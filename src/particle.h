#ifndef _PARTICLE_H
#define _PARTICLE_H
#include <math.h>
#include <stdio.h>

enum ParticleFate  {moving, stopped, escaped, absorbed};

// particle class
class particle
{

public:
  
  double    x[3];         // x,y,z position
  double    D[3];         // direction vector, Dx,Dy,Dz
  int        ind;         // index of the zone in grid where we are
  double       t;         // current time
  double       e;         // total energy in ergs of packet
  double      nu;         // frequency
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

  void print() const
  {
    printf("%10.e | %10.3e %10.3e %10.3e %6d | %10.3e %10.3e %10.3e | %10.3e %10.3e\n",
	   t,x[0],x[1],x[2],ind,D[0],D[1],D[2],e,nu);
  }

};

#endif
