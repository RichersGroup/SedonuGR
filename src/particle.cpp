#include "global_options.h"
#include "particle.h"

particle::particle(){
	x.resize(3,NaN);
	D.resize(3,NaN);
	t = NaN;
	e = NaN;
	nu = NaN;
	s = -1;
	fate = stopped;
}
