#ifndef _EINSTEINHELPER_H
#define _EINSTEINHELPER_H 1

#include "global_options.h"
#include "Particle.h"
#include "Metric.h"
#include "physical_constants.h"
#include "MultiDArray.h"

namespace pc = physical_constants;
using namespace std;

enum TetradRotation {cartesian, spherical};

class EinsteinHelper{
public:
	// essential variables interpolated from grid
	Tuple<double,4> xup;
	Tuple<double,4> kup, kup_tet; // erg
	Tuple<double,4> u; // dimensionless, up index
	Tuple<double,3> v; // cm/s
	double N, s;
	ParticleFate fate;
	double N0;
	Metric g;
	Christoffel christoffel;

	// things with which to do interpolation
	InterpolationCube<NDIMS  > icube_vol; // for metric quantities
	InterpolationCube<NDIMS+1> icube_spec; // for eas

	// intermediate quantities
	Tuple<double,4> e[4]; // [tet(low)][coord(up)]
	double grid_coords[NDIMS+1];
	double absopac, scatopac;
	double ds_com;
	unsigned dir_ind[NDIMS+1]; // spatial, nu_in
	int z_ind, eas_ind;   // direct access indices

	EinsteinHelper(){
		for(unsigned i=0; i<4; i++){
			xup[i] = NaN;
			kup[i] = NaN;
			u[i] = NaN;
			kup_tet[i] = NaN;
			for(unsigned j=0; j<4; j++)
				e[i][j] = NaN;
		}
		N0 = NaN;
		N = NaN;
		s = -MAXLIM;
		fate = moving;
		absopac = NaN;
		scatopac = NaN;
		ds_com = NaN;
		z_ind = -MAXLIM;
		eas_ind = -MAXLIM;
		for(unsigned i=0; i<NDIMS+1; i++)
			dir_ind[i] = MAXLIM;
		for(unsigned i=0; i<3; i++) v[i] = NaN;
	}

	void set_kup_tet(const Tuple<double,4>& kup_tet_in){
		PRINT_ASSERT(Metric::dot_Minkowski<4>(kup_tet_in,kup_tet_in)/(kup_tet_in[3]*kup_tet_in[3]),<,TINY);
		for(unsigned i=0; i<4; i++) kup_tet[i] = kup_tet_in[i];
		tetrad_to_coord(kup_tet,kup);
		g.normalize_null_preserveupt(kup);
		PRINT_ASSERT(g.dot<4>(kup,kup)/(kup[3]*kup[3]),<,TINY);
	}
	void renormalize_kup(){
		g.normalize_null_preserveupt(kup);
		coord_to_tetrad(kup, kup_tet);
		PRINT_ASSERT(kup_tet[3],>,0);
		Metric::normalize_null_Minkowski(kup_tet);
	}

	// return the Lorentz factor W
	static double lorentzFactor(const Metric* g, const Tuple<double,3>& v){
		double result = 1. / sqrt(1. - g->dot<3>(v,v));
		PRINT_ASSERT(result,>=,1);
		return result;
	}

	double nu() const{
		double nu = kup_tet[3] / pc::h;
		PRINT_ASSERT(nu,>,0);
		return nu;
	}

	// get four velocity from three velocity
	void set_fourvel(){
		const Tuple<double,3> vdimless = {v[0]/pc::c, v[1]/pc::c, v[2]/pc::c};
		double W = lorentzFactor(&g,vdimless);
		u[3] = W / (DO_GR ? g.alpha : 1.0);
		PRINT_ASSERT(u[3],>,0);
		PRINT_ASSERT(u[3],<,INFINITY);
		for(unsigned i=0; i<3; i++){
			u[i] = W*vdimless[i];
			if(DO_GR) u[i] -= W/g.alpha * g.betaup[i];
		}
		PRINT_ASSERT(fabs(g.dot<4>(u,u)+1.0),<,TINY);
	}


	// get a Cartesian tetrad basis
	void set_tetrad_basis(TetradRotation rotation){
	  // set the tetrad guesses
	  if(rotation == cartesian){
	    e[0][0] = 1.0;
	    e[0][1] = 0;
	    e[0][2] = 0;
	    e[0][3] = 0;
	    
	    e[1][0] = 0;
	    e[1][1] = 1.0;
	    e[1][2] = 0;
	    e[1][3] = 0;
	    
	    e[2][0] = 0;
	    e[2][1] = 0;
	    e[2][2] = 1.0;
	    e[2][3] = 0;
	  }

	  else if(rotation == spherical){
	    const double rp = sqrt(xup[0]*xup[0] + xup[1]*xup[1]);
	    
	    // pathological case
	    if(rp==0){
	      e[0][0] = 1.0; // theta in x direction
	      e[1][1] = 1.0; // phi in y direction
	      e[2][2] = xup[2]>0 ? 1.0 : -1.0; // radial in z direction
	    }
	    else{
	      // theta vector (multiplying all components by rp*r so no divide by zeros)
	      e[0][0] = xup[0] * xup[2];
	      e[0][1] = xup[1] * xup[2];
	      e[0][2] = -rp * rp;
	      e[0][3] = 0;
	      
	      // phi vector (multiplying through by rp so no divide by zero)
	      e[1][0] = -xup[1];
	      e[1][1] =  xup[0];;
	      e[1][2] = 0;
	      e[1][3] = 0;

	      // radial vector
	      e[2][0] = xup[0];
	      e[2][1] = xup[1];
	      e[2][2] = xup[2];
	      e[2][3] = 0;
	    }
	  }
	  else assert(0);

	  // normalize four-velocity to get timelike vector
	  for(int mu=0; mu<4; mu++) e[3][mu] = u[mu];
	  g.normalize(e[3]);
	  
	  // use x0 as a trial vector
	  g.orthogonalize<4>(e[2],e[3]);
	  g.normalize(e[2]);

	  // use x1 as a trial vector
	  g.orthogonalize<4>(e[1],e[3]);
	  g.orthogonalize<4>(e[1],e[2]);
	  g.normalize(e[1]);
	
	  // use x2 as a trial vector
	  g.orthogonalize<4>(e[0],e[3]);
	  g.orthogonalize<4>(e[0],e[2]);
	  g.orthogonalize<4>(e[0],e[1]);
	  g.normalize(e[0]);
	  
	  // sanity checks
	  PRINT_ASSERT(fabs(g.dot<4>(e[0],e[1])),<,TINY);
	  PRINT_ASSERT(fabs(g.dot<4>(e[0],e[2])),<,TINY);
	  PRINT_ASSERT(fabs(g.dot<4>(e[0],e[3])),<,TINY);
	  PRINT_ASSERT(fabs(g.dot<4>(e[1],e[2])),<,TINY);
	  PRINT_ASSERT(fabs(g.dot<4>(e[1],e[3])),<,TINY);
	  PRINT_ASSERT(fabs(g.dot<4>(e[2],e[3])),<,TINY);
	}

	void coord_to_tetrad(const Tuple<double,4>& kup_coord, Tuple<double,4>& kup_tet) const{
		for(int mu=0; mu<4; mu++) kup_tet[mu] = g.dot<4>(kup_coord,e[mu]);
		kup_tet[3] *= -1.; // k.e = kdown_tet. Must raise index.
	}

	void tetrad_to_coord(const Tuple<double,4>& kup_tet, Tuple<double,4>& kup_coord) const{
		for(int mu=0; mu<4; mu++){
			kup_coord[mu] = 0;
			for(int nu=0; nu<4; nu++)
				kup_coord[mu] += kup_tet[nu] * e[nu][mu];
		}
	}

	template<unsigned n>
	double dot_tetrad(const Tuple<double,4>& x1, const Tuple<double,4>& x2){
		double result = 0;
		for(unsigned i=0; i<3; i++) result += x1[n]*x2[n];
		if(n==4) result -= x1[2]*x2[3];
		return result;
	}

	void scale_p_frequency(const double scale){
		kup *= scale;
	}

	void get_Particle(ParticleList& pout, const unsigned list_index) const{
		pout.N[list_index] = N;
		pout.s[list_index] = s;
		pout.fate[list_index] = fate;
		for(unsigned i=0; i<4; i++){
			pout.xup[i][list_index] = xup[i];
			pout.kup[i][list_index] = kup[i];
		}
	}
	void set_Particle(const ParticleList& pin, const unsigned list_index){
		N = pin.N[list_index];
		s = pin.s[list_index];
		fate = pin.fate[list_index];
		for(unsigned i=0; i<4; i++){
			xup[i] = pin.xup[i][list_index];
			kup[i] = pin.kup[i][list_index];
		}
	}
};

#endif
