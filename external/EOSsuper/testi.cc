// -*-C++-*-

#include "vecmathlib.h"
#include <iostream>
#include <vector>
#include <cstdlib>
#include <cstdio>

using namespace std;
using namespace vecmathlib;

int main(void) {
  // Choose an "interesting" vector type
  typedef realvec<double,2> realvec_t;

  
  clock_t start1, end1;
  clock_t start2, end2;
  double elapsed1,elapsed2;
  start1 = clock();
  // set up arrays for testing EOS performance
  double* __restrict__ rhos; 
  double* __restrict__ temps;
  double* __restrict__ yets;
  double* __restrict__ xres;
  const int npoints = 40000000;
  const double lrhomin = 6.0;
  const double lrhomax = 15.0;
  const double ltempmin = 0.05;
  const double ltempmax = 1.7;
  const double yemin = 0.05;
  const double yemax = 0.52;

  rhos = (double*)malloc(sizeof(double)*npoints);
  temps = (double*)malloc(sizeof(double)*npoints);
  yets = (double*)malloc(sizeof(double)*npoints);
  xres = (double*)malloc(sizeof(double)*npoints);
  double* __restrict__ xres2 = (double*)malloc(sizeof(double)*npoints);
  int* __restrict__ keyerr = (int*)malloc(sizeof(int)*npoints);
  int anyerr;

  const int seed = 1234567;
  unsigned int rand_state = (unsigned int)seed;
  for(int i=0;i<npoints;i++) {

    const double xrho = pow(lrhomin + (lrhomax-lrhomin) * 1.0 *
			    rand_r(&rand_state)/((double)RAND_MAX),10);
    const double xtemp = pow(ltempmin + (ltempmax-ltempmin) * 1.0 *
			     rand_r(&rand_state)/((double)RAND_MAX),10);
    const double xye = yemin + (yemax-yemin) * 1.0 *
      rand_r(&rand_state)/((double)RAND_MAX);

    rhos[i] = xrho;
    temps[i] = xtemp;
    yets[i] = xye;
    xres[i] = xrho;
    xres2[i] = xye;
  }
  end1 = clock();
  elapsed1 = ((double) (end1 - start1)) / CLOCKS_PER_SEC;


#if 1
  typedef realvec<double,2>  double2;
  vector<double2> xvres(npoints/2);
  for(int i=0;i<npoints/2;++i) xvres[i] = double2(&xres[2*i]);
  vector<double2> xvres2(npoints/2);
  for(int i=0;i<npoints/2;++i) xvres2[i] = double2(&xres2[2*i]);
#endif

  
  start2 = clock();
#if 0
  for(int i=0;i<npoints;i++) {
    //xres[i] = pow(10.0,xres[i]);
    //xres2[i] = pow(10.0,xres2[i]);
    xres[i] = exp(M_LN10 * xres[i]);
    xres2[i] = exp(M_LN10 * xres2[i]);
  }
#endif

#if 1
  for (ptrdiff_t j=0; j<npoints/2; ++j) {
    xvres[j] = exp(double2(M_LN10)*xvres[j]);
    xvres2[j] = exp(double2(M_LN10)*xvres2[j]);
  }
#endif

  end2 = clock();
  elapsed2 = ((double) (end2 - start2)) / CLOCKS_PER_SEC;

  fprintf(stderr,"%8.5f %8.5f\n",
	  elapsed1,elapsed2);



  return 0;
}
