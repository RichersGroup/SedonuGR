#ifndef _METRIC_H
#define _METRIC_H 1

#include "global_options.h"
#include "gsl/gsl_linalg.h"

const size_t ixx=0,iyy=1,izz=2,ixy=3,ixz=4,iyz=5,itt=6,ixt=7,iyt=8,izt=9;

class ThreeMetric{
// both indices are down
public:
	Tuple<double,6> data;

	static int index(const size_t i, const size_t j){
		PRINT_ASSERT(i,<,3);
		PRINT_ASSERT(j,<,3);
		switch( (j+1)*(i+1) ){
		case 1:
			return ixx;
		case 2:
			return ixy;
		case 3:
			return ixz;
		case 4:
			return iyy;
		case 6:
			return iyz;
		case 9:
			return izz;
		default:
			return -1;
		}
	}

 ThreeMetric() : data(NaN) {}

	double det() const{
		double result = 0;
		result += data[ixx] * (data[iyy]*data[izz] - data[iyz]*data[iyz]);
		result -= data[ixy] * (data[ixy]*data[izz] - data[iyz]*data[ixz]);
		result += data[ixz] * (data[ixy]*data[iyz] - data[iyy]*data[ixz]);
		PRINT_ASSERT(result,>,0);
		return result;
	}
	Tuple<double,3> lower(const Tuple<double,3>& in) const{
	  Tuple<double,3> out;
	  out[0] = in[0]*data[ixx] + in[1]*data[ixy] + in[2]*data[ixz];
	  out[1] = in[0]*data[ixy] + in[1]*data[iyy] + in[2]*data[iyz];
	  out[2] = in[0]*data[ixz] + in[1]*data[iyz] + in[2]*data[izz];
	  return out;
	}
	ThreeMetric inverse() const{
		gsl_matrix* g = gsl_matrix_alloc(3,3);
		for(size_t i=0; i<3; i++)
		  for(size_t j=0; j<3; j++)
		    gsl_matrix_set(g,i,j, data[index(i,j)]);

		// invert and store
		gsl_linalg_cholesky_decomp(g);
		gsl_linalg_cholesky_invert(g);
		ThreeMetric output;
		output.data[ixx] = gsl_matrix_get(g, 0, 0);
		output.data[ixy] = gsl_matrix_get(g, 0, 1);
		output.data[ixz] = gsl_matrix_get(g, 0, 2);
		output.data[iyy] = gsl_matrix_get(g, 1, 1);
		output.data[iyz] = gsl_matrix_get(g, 1, 2);
		output.data[izz] = gsl_matrix_get(g, 2, 2);

		// free the memory
		gsl_matrix_free(g);

		return output;
	}

	inline double get(const size_t i, const size_t j) const{
		return data[index(i,j)];
	}
};

//========//
// METRIC //
//========//
class Metric{
public:
	Tuple<double,3> betalow, betaup;
	double gtt, alpha;
	ThreeMetric gammalow, gammaup;

 Metric() : betalow(NaN), betaup(NaN), gtt(NaN), alpha(NaN){}

	static int index(const size_t i, const size_t j){
		PRINT_ASSERT(i,<,4);
		PRINT_ASSERT(j,<,4);
		switch( (j+1)*(i+1) ){
		case 1:
			return ixx;
		case 2:
			return ixy;
		case 3:
			return ixz;
		case 4:
			return (i+j==2 ? iyy : ixt);
		case 6:
			return iyz;
		case 9:
			return izz;
		case 8:
			return iyt;
		case 12:
			return izt;
		case 16:
			return itt;
		default:
			return -1;
		}
	}

	// fill in values for gtt and betalow
	// assumes alpha, betaup, and gammalow have been set.
	void update(){
	  assert(DO_GR);
	  betalow = gammalow.lower(betaup);
	  gtt = DO_GR ? -alpha*alpha + contract<3>(betaup, betalow) : -1.0;
	  gammaup = gammalow.inverse();
	}

	double get(const size_t i, const size_t j) const{
		if(i==3 and j==3) return gtt;
		else if(i==3) return betalow[j];
		else if(j==3) return betalow[i];
		else return gammalow.get(i,j);
	}
	double get_inverse(const size_t i, const size_t j) const{
		PRINT_ASSERT(gammaup.get(0,0),==,gammaup.get(0,0));
		if(i==3 and j==3) return -1./(alpha*alpha);
		else if(i==3) return betaup[j]/(alpha*alpha);
		else if(j==3) return betaup[i]/(alpha*alpha);
		else return gammaup.get(i,j) - betaup[i]*betaup[j]/(alpha*alpha);
	}

	template<size_t n, size_t n1>
	Tuple<double,n> lower(const Tuple<double,n1>& xup) const{
		PRINT_ASSERT(n,<=,n1);
		Tuple<double,n> xdown;
		if(DO_GR){
			for(size_t i=0; i<n; i++){
				xdown[i] = 0;
				for(size_t j=0; j<n; j++)
					xdown[i] += xup[j] * get(i,j);
			}
		}
		else{
			for(size_t i=0; i<3; i++) xdown[i] = xup[i];
			if(n==4) xdown[3] = -xup[3];
		}
		return xdown;
	}

	template<size_t n>
	Tuple<double,n> raise(const Tuple<double,n>& xdown) const{
	        Tuple<double,n> xup;
		if(DO_GR){
			for(size_t i=0; i<n; i++){
				xup[i] = 0;
				for(size_t j=0; j<n; j++)
					xup[i] += xdown[j] * get_inverse(i,j);
			}
		}
		else{
			for(size_t i=0; i<3; i++) xup[i] = xdown[i];
			if(n==4) xup[3] = -xdown[3];
		}
		return xup;
	}

	template<size_t n, size_t n1, size_t n2>
	static double contract(const Tuple<double,n1>& xup, const Tuple<double,n2>& xdown){
		PRINT_ASSERT(n,<=,n1);
		PRINT_ASSERT(n,<=,n2);
		double result = 0;
		for(size_t i=0; i<n; i++) result += xup[i]*xdown[i];
		return result;
	}

	// dot product
	template<size_t n, size_t n1, size_t n2>
	double dot(const Tuple<double,n1>& x1up, const Tuple<double,n2>& x2up) const{
		PRINT_ASSERT(n,<=,n1);
		PRINT_ASSERT(n,<=,n2);
		if(DO_GR){
			Tuple<double,n> x2low = lower<n>(x2up);
			double result = contract<n>(x1up, x2low);
			return result;
		}
		else return dot_Minkowski<n>(x1up, x2up);
	}

	// dot the normal observer's four-velocity with a four vector
	double ndot(const Tuple<double,4>& x) const{
		return -(DO_GR ? alpha : 1.0) * x[3];
	}

	// normalize a four vector to have a norm of +/-1
	void normalize(Tuple<double,4>& x) const{
		double invnorm = sqrt(fabs(1./dot<4>(x,x)));
		PRINT_ASSERT(invnorm,>,0);
		for(size_t i=0; i<4; i++) x[i] *= invnorm;
	}

	// make a vector null
	void normalize_null_preservedownt(Tuple<double,4>& x) const{
		PRINT_ASSERT(x[3],>=,0);
		double result = NaN;
		if(DO_GR){
			Tuple<double,4> xlow = lower<4>(x);
			double C = get_inverse(3,3) * xlow[3]*xlow[3];
			double B=0, A=0;
			for(size_t i=0; i<3; i++){
				B += 2.*get_inverse(i,3) * xlow[3]*xlow[i];
				for(size_t j=0; j<3; j++)
					A += get_inverse(i,j) * xlow[i]*xlow[j];
			}

			B /= A;
			C /= A;
			PRINT_ASSERT(B*B - 4.*C,>=,0);
			result = 0.5 * (-B + sqrt(B*B - 4.*C));
			PRINT_ASSERT(result,>,TINY);

			for(size_t i=0; i<3; i++) xlow[i] *= result;
			x = raise(xlow);
			PRINT_ASSERT(dot<4>(x,x)/(x[3]*x[3]),<,TINY);
			PRINT_ASSERT(x[3],>=,0);
		}
		else{
			result = x[3] / sqrt(dot_Minkowski<3>(x,x));
			for(size_t i=0; i<3; i++) x[i] *= result;
		}
		PRINT_ASSERT(abs(dot<4>(x,x))/(x[3]*x[3]),<,TINY);
		PRINT_ASSERT(x[3],>=,0);
	}

	void normalize_null_preserveupt(Tuple<double,4>& x) const{
		double result = NaN;
		if(DO_GR){
			double A = dot<3>(x,x);
			double C = gtt * x[3]*x[3] / A;
			double B = 2.*x[3] * contract<3>(betalow,x) / A;
			PRINT_ASSERT(B*B - 4.*C,>=,0);
			result = 0.5 * (-B + sqrt(B*B - 4.*C));

			for(size_t i=0; i<3; i++) x[i] *= result;
			PRINT_ASSERT(dot<4>(x,x)/(x[3]*x[3]),<,TINY);
		}
		else{
			result = x[3] / sqrt(dot_Minkowski<3>(x,x));
			for(size_t i=0; i<3; i++) x[i] *= result;
		}
		PRINT_ASSERT(abs(dot<4>(x,x))/(x[3]*x[3]),<,TINY);
	}

	void normalize_null_changeupt(Tuple<double,4>& x) const{
		double result = NaN;
		if(DO_GR){
			const double invA = 1./gtt;
			const double B = 2.*contract<3>(betalow,x) * invA;
			const double C = dot<3>(x,x) * invA;
			PRINT_ASSERT(B*B - 4.*C,>=,0);
			result = 0.5 * (-B + sqrt(B*B - 4.*C));
		}
		else result = sqrt(dot_Minkowski<3>(x,x));
		x[3] = result;
		PRINT_ASSERT(x[3],>=,0);
		PRINT_ASSERT(abs(dot<4>(x,x))/(x[3]*x[3]),<,TINY);
	}

	// make four vector v orthogonal to four vector v
	template<size_t n>
	void orthogonalize(Tuple<double,n>& v, const Tuple<double,n>& e) const{
	        double edote = dot<n>(e,e);
		PRINT_ASSERT(edote,!=,0);
		double projection = dot<n>(v,e) / edote;
		for(size_t mu=0; mu<n; mu++) v[mu] -= projection * e[mu];
	}

	// vector operations
	template<size_t s, size_t n1, size_t n2>
	static double dot_Minkowski(const Tuple<double,n1>& a, const Tuple<double,n2>& b){
		PRINT_ASSERT(s,<=,n1);
		PRINT_ASSERT(s,<=,n2);
		double product = 0;
		for(size_t i=0; i<3; i++) product += a[i]*b[i];
		if(s==4) product -= a[3]*b[3];
		return product;
	}

	// normalize a vector
	template<size_t s>
	static void normalize_Minkowski(Tuple<double,s>& a){
		double inv_magnitude = 1./sqrt(fabs( dot_Minkowski<s>(a,a) ));
		PRINT_ASSERT(inv_magnitude,<,INFINITY);
		for(size_t i=0; i<s; i++) a[i] *= inv_magnitude;
	}

	static void normalize_null_Minkowski(Tuple<double,4>& a){
		double spatial_norm = dot_Minkowski<3>(a,a);
		a[3] = sqrt(spatial_norm);
	}

};



class Christoffel{
public:
	// first index is up, others are down
	// 0-9 = 0 first index, 10-19 = 1 first index, 20-29 = 2 first index, 20-23 - third index
	// 3 first index not included
	Tuple<double, 40> data;

	inline static size_t index(const size_t a, const size_t i, const size_t j){
		PRINT_ASSERT(a,<=,3);
		return 10*a + Metric::index(i,j);
	}

	Christoffel() : data(NaN) {}

	Tuple<double,4> contract2(const Tuple<double,4>& kup) const{
		Tuple<double,4> result = 0;
		for(size_t a=0; a<4; a++){
			for(size_t i=0; i<4; i++)
				for(size_t j=0; j<4; j++)
					result[a] += data[index(a,i,j)] * kup[i]*kup[j];
		}
		return result;
	}
};
#endif
