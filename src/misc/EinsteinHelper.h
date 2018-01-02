#ifndef _EINSTEINHELPER_H
#define _EINSTEINHELPER_H 1

#include "global_options.h"
#include "Particle.h"
#include "Metric.h"
#include "physical_constants.h"

namespace pc = physical_constants;
using namespace std;

class EinsteinHelper{
public:
  // essential variables interpolated from grid
  Particle p;
  Metric g;
  double u[4]; // dimensionless

  // intermediate quantities
  double e[4][4];
  double kup_tet[4];
  double absopac, scatopac;
  double ds_com;
  vector<unsigned> dir_ind;

  // fill in values for g.{gtt, betalow}, u, e, kup_tet
  // assumes g.{alpha, betaup, gammalow}, p.kup are set
  void update(const double v[3]){ // use v as input since we don't store it
    g.update();
    set_fourvel(v);
    set_tetrad_basis(u);
    coord_to_tetrad(p.kup, kup_tet);
    g.normalize_null(p.kup);
  }
  
  // return the Lorentz factor W
  double lorentzFactor(const double v[3]) const{
    return 1. / sqrt(1. - g.dot<3>(v,v));
  }
  
  double nu() const{
	  return -g.dot<4>(p.kup,u) * pc::c / (2.0*pc::pi);
  }

  // get four velocity from three velocity
  void set_fourvel(const double v[3]){
    double W = lorentzFactor(v);
    u[3] = W/g.alpha;
    for(unsigned i=0; i<3; i++)
      u[i] = v[3]/pc::c * (g.alpha*v[i]/pc::c - g.betaup[i]);
  }


  // get a Cartesian tetrad basis
  void set_tetrad_basis(const double u[4]){
    // normalize four-velocity to get timelike vector
    for(int mu=0; mu<4; mu++) e[3][mu] = u[mu];
    g.normalize(e[3]);
    
    // use x0 as a trial vector
    e[0][0] = 1.0;
    e[0][1] = 0;
    e[0][2] = 0;
    e[0][3] = 0;
    g.orthogonalize<4>(e[0],e[3]);
    g.normalize(e[0]);
    
    // use x1 as a trial vector
    e[1][0] = 0;
    e[1][1] = 1.0;
    e[1][2] = 0;
    e[1][3] = 0;
    g.orthogonalize<4>(e[1],e[3]);
    g.orthogonalize<4>(e[1],e[0]);
    g.normalize(e[1]);
    
    // use x2 as a trial vector
    e[2][0] = 0;
    e[2][1] = 0;
    e[2][2] = 1.0;
    e[2][3] = 0;
    g.orthogonalize<4>(e[2],e[3]);
    g.orthogonalize<4>(e[2],e[0]);
    g.orthogonalize<4>(e[2],e[1]);
    g.normalize(e[2]);
    
    // sanity checks
    PRINT_ASSERT(abs(g.dot<4>(e[0],e[1])),<,TINY);
    PRINT_ASSERT(abs(g.dot<4>(e[0],e[2])),<,TINY);
    PRINT_ASSERT(abs(g.dot<4>(e[0],e[3])),<,TINY);
    PRINT_ASSERT(abs(g.dot<4>(e[1],e[2])),<,TINY);
    PRINT_ASSERT(abs(g.dot<4>(e[1],e[3])),<,TINY);
    PRINT_ASSERT(abs(g.dot<4>(e[2],e[3])),<,TINY);
  }
  
  void coord_to_tetrad(const double kup_coord[4], double kup_tet[4]) const{
    for(int mu=0; mu<4; mu++) kup_tet[mu] = g.dot<4>(kup_coord,e[mu]);
  }

  void tetrad_to_coord(const double kup_tet[4], double kup_coord[4]) const{
    for(int mu=0; mu<4; mu++){
      kup_coord[mu] = 0;
      for(int nu=0; nu<4; nu++)
	kup_coord[mu] += kup_tet[nu] * e[nu][mu];
    }
  }

  template<unsigned n>
    double dot_tetrad(const double x1[n], const double x2[n]){
    double result = 0;
    for(unsigned i=0; i<3; i++) result += x1[n]*x2[n];
    if(n==4) result -= x1[2]*x2[3];
    return result;
  }

  void scale_p_frequency(const double scale){
	  for(unsigned i=0; i<4; i++) p.kup[i] *= scale;
  }

  double ds_lab(const double ds_com_in){
	  return ds_com_in * g.ndot(p.kup)/g.dot<4>(u,p.kup);
  }
};

#endif
