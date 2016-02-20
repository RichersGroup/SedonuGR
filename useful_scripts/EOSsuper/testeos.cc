#include <stdio.h>
#include <math.h>
#include <time.h>
#include "nuc_eos.hh"
#include "helpers.hh"
//#include "vecmathlib.h"
#include <iostream>
#include <vector>
#include <cstdlib>
#define restrict __restrict__


using namespace std;
//using namespace vecmathlib;

namespace nuc_eos {

void testeos(void) {

  using namespace nuc_eos;

  // Choose a vector size
#if defined VECMATHLIB_HAVE_VEC_DOUBLE_4
  //  typedef realvec<double,4> realvec_t;
#elif defined VECMATHLIB_HAVE_VEC_DOUBLE_2
  //typedef realvec<double,2> realvec_t;
#else
  //typedef realpseudovec<double,1> realvec_t;
#endif

  clock_t start1, end1;
  clock_t start2, end2;
  double elapsed1,elapsed2;
  start1 = clock();
  // set up arrays for testing EOS performance
  double* __restrict__ rhos; 
  double* __restrict__ temps;
  double* __restrict__ yets;
  double* __restrict__ xres;
  const int npoints = 10000000;
  const double lrhomin = log(1.0e6);
  const double lrhomax = log(1.0e15);
  const double ltempmin = log(0.03);
  const double ltempmax = log(80);
  const double yemin = 0.05;
  const double yemax = 0.52;

  rhos = (double*)malloc(sizeof(double)*npoints);
  temps = (double*)malloc(sizeof(double)*npoints);
  yets = (double*)malloc(sizeof(double)*npoints);
  xres = (double*)malloc(sizeof(double)*npoints);
  double* __restrict__ xres2 = (double*)malloc(sizeof(double)*npoints);
  double* __restrict__ xres3 = (double*)malloc(sizeof(double)*npoints);
  int* __restrict__ keyerr = (int*)malloc(sizeof(int)*npoints);
  int anyerr;

  const int seed = 1234567;
  unsigned int rand_state = (unsigned int)seed;
  for(int i=0;i<npoints;i++) {

    const double xrho = exp(lrhomin + (lrhomax-lrhomin) * 1.0 *
			    rand_r(&rand_state)/((double)RAND_MAX));
    const double xtemp = exp(ltempmin + (ltempmax-ltempmin) * 1.0 *
    			     rand_r(&rand_state)/((double)RAND_MAX));
    const double xye = yemin + (yemax-yemin) * 1.0 *
      rand_r(&rand_state)/((double)RAND_MAX);

    rhos[i] = xrho*RHOGF;
    temps[i] = xtemp;
    yets[i] = xye;
    xres[i] = 0.0e0;
    xres2[i] = 0.0e0;
  }
  end1 = clock();
  elapsed1 = ((double) (end1 - start1)) / CLOCKS_PER_SEC;


  nuc_eos_m_kt1_press_eps(&npoints,rhos,temps,yets,xres,xres2,
			  keyerr,&anyerr);

#if 0
  typedef realvec<double,2>  double2;
  vector<double2> xvres(npoints/2);
  for(int i=0;i<npoints;i+=2) xvres[i/2] = double2(&xres[i]);
  vector<double2> xvres2(npoints/2);
  for(int i=0;i<npoints;i+=2) xvres2[i/2] = double2(&xres2[i]);
#endif

#if 0
  int j=116856;
  fprintf(stderr,"%15.6E %22.14E %15.6E %15.6E %15.6E\n",
	  rhos[j],temps[j],yets[j],xres2[j],xres[j]);
#endif    

  for(int i=0;i<npoints;i++) 
  { 
    double rf = 1.0e-1 * rand_r(&rand_state)/((double)RAND_MAX);
    //    xres[i] += xres[i]*rf;
    //    fprintf(stderr,"%d, %15.6E\n",i,xres[i]);
    temps[i] += temps[i]*rf; 
    temps[i] = MIN(100.0,temps[i]);
    //temps[i] = 0.014e0;
  }

#if 0
  fprintf(stderr,"%15.6E %22.14E %15.6E %15.6E %15.6E\n",
	  rhos[j],temps[j],yets[j],xres2[j],xres[j]);
#endif

  start2 = clock();
  double prec = 1.0e-12;
  nuc_eos_m_kt0_press(&npoints,rhos,temps,yets,
		      xres,xres2,&prec,keyerr,&anyerr);
  if(anyerr) {
    fprintf(stderr,"Caught exception. Something bad happening!\n");
    int i=0;
    while(i<npoints) {
      if(keyerr[i]) {
	fprintf(stderr,
		"i: %d   keyerr: %d  rho: %15.6E  t: %15.6E  ye: %15.6E eps: %15.6E\n",
		i,keyerr[i],
		rhos[i],temps[i],yets[i],xres[i]);
      }
      i++;
    }
  }
#if 0
  fprintf(stderr,"%15.6E %22.14E %15.6E %15.6E %15.6E\n",
	  rhos[j],temps[j],yets[j],xres2[j],xres[j]);
#endif

  end2 = clock();
  elapsed2 = ((double) (end2 - start2)) / CLOCKS_PER_SEC;

  fprintf(stderr,"%8.5f %8.5f\n",
	  elapsed1,elapsed2);

  return;
}


} // namespace


#if 0
  for(int i=0;i<npoints;i++) {
    xres[i] = exp(xres[i]);
    xres2[i] = exp(xres2[i]);
  }
#endif

#if 0
  for (ptrdiff_t j=0; j<npoints/2; ++j) {
    //    xvres[j] = exp(M_LN10*xres[j]);
    //xvres2[j] = exp(M_LN10*xres2[j]);
  }
#endif
