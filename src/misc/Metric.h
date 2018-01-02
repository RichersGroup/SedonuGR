#ifndef _METRIC_H
#define _METRIC_H 1

#include "global_options.h"

class ThreeMetric{
public:
  double xx, yy, zz, xy, xz, yz;
};

//========//
// METRIC //
//========//
class Metric{
 public:
  double gtt, alpha, betalow[3], betaup[3];
  ThreeMetric gammalow;

  // fill in values for gtt and betalow
  // assumes alpha, betaup, and gammalow have been set.
  void update(){
    lower<4>(betaup, betalow);
    gtt = -alpha*alpha + contract<4>(betaup, betalow);
  }

  template<unsigned n>
  void lower(const double xup[], double xdown[]) const{
    xdown[0] = gammalow.xx*xup[0] + gammalow.xy*xup[1] + gammalow.xz*xup[2];
    xdown[1] = gammalow.xy*xup[0] + gammalow.yy*xup[1] + gammalow.yz*xup[2];
    xdown[2] = gammalow.xz*xup[0] + gammalow.yz*xup[1] + gammalow.zz*xup[2];

    if(n==4){
      xdown[3]  = xup[3] * gtt;
      for(unsigned i=0; i<3; i++){
	xdown[3] += xup[i] * betalow[i];
	xdown[i] += xup[3] * betalow[i];
      }
    }
  }

  template<unsigned n>
    double contract(const double xup[n], const double xdown[n]) const{
    double result = 0;
    for(int i=0; i<n; i++) result += xup[n]*xdown[n];
    return result;
  }
        
  // dot product
  template<unsigned n>
  double dot(const double x1up[n], const double x2up[n]) const{
    double x2low[n];
    lower<n>(x2up, x2low);
    double result = contract<n>(x1up, x2low);
    return result;
  }

  // dot the normal observer's four-velocity with a four vector
  double ndot(const double x[4]) const{
    return -alpha * x[3];
  }

  // normalize a four vector to have a norm of +/-1
  void normalize(double x[4]) const{
    double invnorm = abs(1./dot<4>(x,x));
    for(int i=0; i<4; i++) x[i] *= invnorm;
  }

  // make a vector null
  void normalize_null(double x[4]) const{
    const double invA = 1./gtt;
    const double B = (x[0]*betalow[0] + x[1]*betalow[1] + x[2]*betalow[2]) * invA;
    const double C = dot<3>(x,x) * invA;
    const double result = 0.5 * (B + sqrt(B*B - 4.*C));
    PRINT_ASSERT(result,>,0);
    x[3] = result;
  }

  // make four vector v orthogonal to four vector v
  template<unsigned n>
  void orthogonalize(double v[n], const double e[n]) const{
    double projection = dot<n>(v,e);
    for(int mu=0; mu<n; mu++) v[mu] -= projection * e[mu];
  }  
};

#endif
