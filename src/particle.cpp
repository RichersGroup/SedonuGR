#include "particle.h"
#include <limits>
#define NaN std::numeric_limits<double>::quiet_NaN()

particle::particle(){
	x.resize(3,NaN);
	D.resize(3,NaN);
	t = NaN;
	e = NaN;
	nu = NaN;
	s = NaN;
	fate = stopped;
}

void particle::normalize_direction(){
    double magD = sqrt(D[0]*D[0] + D[1]*D[1] + D[2]*D[2]);
    D[0] /= magD;
    D[1] /= magD;
    D[2] /= magD;
}
