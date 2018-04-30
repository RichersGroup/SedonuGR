#ifndef _METRIC_H
#define _METRIC_H 1

#include "global_options.h"
#include "gsl/gsl_linalg.h"

const unsigned ixx=0,iyy=1,izz=2,ixy=3,ixz=4,iyz=5,itt=6,ixt=7,iyt=8,izt=9;

class ThreeMetric{
// both indices are down
public:
	Tuple<double,6> data;

	static int index(const unsigned i, const unsigned j){
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
			return (i+j==2 ? iyy : ixt);
		case 6:
			return iyz;
		case 9:
			return izz;
		default:
			return -1;
		}
	}

	ThreeMetric(){
		data = NaN;
	}

	double det() const{
		double result = 0;
		result += data[ixx] * (data[iyy]*data[izz] - data[iyz]*data[iyz]);
		result -= data[ixy] * (data[ixy]*data[izz] - data[iyz]*data[ixz]);
		result += data[ixz] * (data[ixy]*data[iyz] - data[iyy]*data[ixz]);
		return result;
	}
	void lower(const Tuple<double,3>& in, Tuple<double,3>& out) const{
		for(unsigned i=0; i<3; i++){
			out[i] = 0;
			for(unsigned j=0; j<3; j++){
				out[i] += get(i,j)*in[j];
			}
		}
	}
	ThreeMetric inverse() const{
		gsl_matrix* g = gsl_matrix_alloc(3,3);
		for(unsigned i=0; i<3; i++)
			for(unsigned j=0; j<3; j++)
				gsl_matrix_set(g,i,j, data[index(i,j)]);

		// get LU decomposition
		int s;
		gsl_permutation *p = gsl_permutation_alloc(3);
		gsl_linalg_LU_decomp(g, p, &s);

		// invert and store
		gsl_matrix* ginv = gsl_matrix_alloc(3,3);
		gsl_linalg_LU_invert(g, p, ginv);
		ThreeMetric output;
		output.data[ixx] = gsl_matrix_get(ginv, 0, 0);
		output.data[ixy] = gsl_matrix_get(ginv, 0, 1);
		output.data[ixz] = gsl_matrix_get(ginv, 0, 2);
		output.data[iyy] = gsl_matrix_get(ginv, 1, 1);
		output.data[iyz] = gsl_matrix_get(ginv, 1, 2);
		output.data[izz] = gsl_matrix_get(ginv, 2, 2);

		// free the memory
		gsl_permutation_free(p);
		gsl_matrix_free(g);
		gsl_matrix_free(ginv);

		return output;
	}

	double get(const unsigned i, const unsigned j) const{
		return data[index(i,j)];
	}
};

//========//
// METRIC //
//========//
class Metric{
public:
	double gtt, betalow[3];
	double alpha, betaup[3];
	ThreeMetric gammalow, gammaup;

	Metric(){
		gtt=NaN;
		alpha = NaN;
		for(unsigned i=0; i<3; i++){
			betaup[i] = NaN;
			betalow[i] = NaN;
		}
	}

	static int index(const unsigned i, const unsigned j){
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
		if(DO_GR){
			lower<3>(betaup, betalow);
			gtt = DO_GR ? -alpha*alpha + contract<3>(betaup, betalow) : -1.0;
		}
	}
	void set_inverse(){
		if(DO_GR) gammaup = gammalow.inverse();
	}

	double get(const unsigned i, const unsigned j) const{
		if(i==3 and j==3) return gtt;
		else if(i==3) return betalow[j];
		else if(j==3) return betalow[i];
		else return gammalow.get(i,j);
	}
	double get_inverse(const unsigned i, const unsigned j) const{
		PRINT_ASSERT(gammaup.get(0,0),==,gammaup.get(0,0));
		if(i==3 and j==3) return -1./(alpha*alpha);
		else if(i==3) return betaup[j]/(alpha*alpha);
		else if(j==3) return betaup[i]/(alpha*alpha);
		else return gammaup.get(i,j) - betaup[i]*betaup[j]/(alpha*alpha);
	}

	template<unsigned n>
	void lower(const double xup[], double xdown[]) const{
		if(DO_GR){
			for(unsigned i=0; i<n; i++){
				xdown[i] = 0;
				for(unsigned j=0; j<n; j++)
					xdown[i] += xup[j] * get(i,j);
			}
		}
		else{
			for(unsigned i=0; i<3; i++) xdown[i] = xup[i];
			if(n==4) xdown[3] = -xup[3];
		}
	}

	template<unsigned n>
	void raise(const double xdown[], double xup[]) const{
		if(DO_GR){
			for(unsigned i=0; i<n; i++){
				xup[i] = 0;
				for(unsigned j=0; j<n; j++)
					xup[i] += xdown[j] * get_inverse(i,j);
			}
		}
		else{
			for(unsigned i=0; i<3; i++) xup[i] = xdown[i];
			if(n==4) xup[3] = -xdown[3];
		}
	}

	template<unsigned n>
	static double contract(const double xup[n], const double xdown[n]){
		double result = 0;
		for(unsigned i=0; i<n; i++) result += xup[i]*xdown[i];
		return result;
	}

	// dot product
	template<unsigned n>
	double dot(const double x1up[n], const double x2up[n]) const{
		if(DO_GR){
			double x2low[n];
			lower<n>(x2up, x2low);
			double result = contract<n>(x1up, x2low);
			return result;
		}
		else return dot_Minkowski<n>(x1up, x2up);
	}

	// dot the normal observer's four-velocity with a four vector
	double ndot(const double x[4]) const{
		return -alpha * x[3];
	}

	// normalize a four vector to have a norm of +/-1
	void normalize(double x[4]) const{
		double invnorm = sqrt(fabs(1./dot<4>(x,x)));
		for(unsigned i=0; i<4; i++) x[i] *= invnorm;
	}

	// make a vector null
	void normalize_null_preservedownt(double x[4]) const{
		PRINT_ASSERT(x[3],>=,0);
		double result = NaN;
		if(DO_GR){
			double xlow[4];
			lower<4>(x,xlow);
			double C = get_inverse(3,3) * xlow[3]*xlow[3];
			double B=0, A=0;
			for(unsigned i=0; i<3; i++){
				B += 2.*get_inverse(i,3) * xlow[3]*xlow[i];
				for(unsigned j=0; j<3; j++)
					A += get_inverse(i,j) * xlow[i]*xlow[j];
			}

			B /= A;
			C /= A;
			PRINT_ASSERT(B*B - 4.*C,>=,0);
			result = 0.5 * (-B + sqrt(B*B - 4.*C));
			PRINT_ASSERT(result,>,TINY);

			for(unsigned i=0; i<3; i++) xlow[i] *= result;
			raise<4>(xlow, x);
			PRINT_ASSERT(dot<4>(x,x)/(x[3]*x[3]),<,TINY);
			PRINT_ASSERT(x[3],>=,0);
		}
		else{
			result = x[3] / sqrt(dot_Minkowski<3>(x,x));
			for(unsigned i=0; i<3; i++) x[i] *= result;
		}
		PRINT_ASSERT(abs(dot<4>(x,x))/(x[3]*x[3]),<,TINY);
		PRINT_ASSERT(x[3],>=,0);
	}

	void normalize_null_preserveupt(double x[4]) const{
		double result = NaN;
		if(DO_GR){
			double A = dot<3>(x,x);
			double C = gtt * x[3]*x[3] / A;
			double B = 2.*x[3] * contract<3>(betalow,x) / A;
			PRINT_ASSERT(B*B - 4.*C,>=,0);
			result = 0.5 * (-B + sqrt(B*B - 4.*C));

			for(unsigned i=0; i<3; i++) x[i] *= result;
			PRINT_ASSERT(dot<4>(x,x)/(x[3]*x[3]),<,TINY);
		}
		else{
			result = x[3] / sqrt(dot_Minkowski<3>(x,x));
			for(unsigned i=0; i<3; i++) x[i] *= result;
		}
		PRINT_ASSERT(abs(dot<4>(x,x))/(x[3]*x[3]),<,TINY);
	}

	void normalize_null_changeupt(double x[4]) const{
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
	template<unsigned n>
	void orthogonalize(double v[n], const double e[n]) const{
		double projection = dot<n>(v,e) / dot<n>(e,e);
		for(unsigned mu=0; mu<n; mu++) v[mu] -= projection * e[mu];
	}

	// vector operations
	template<unsigned s>
	static double dot_Minkowski(const double a[], const double b[]){
		double product = 0;
		for(unsigned i=0; i<3; i++) product += a[i]*b[i];
		if(s==4) product -= a[3]*b[3];
		return product;
	}

	// normalize a vector
	template<unsigned s>
	static void normalize_Minkowski(double a[]){
		double inv_magnitude = 1./sqrt(fabs( dot_Minkowski<s>(a,a) ));
		PRINT_ASSERT(inv_magnitude,<,INFINITY);
		for(unsigned i=0; i<s; i++) a[i] *= inv_magnitude;
	}

	static void normalize_null_Minkowski(double a[4]){
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

	static unsigned index(const unsigned a, const unsigned i, const unsigned j){
		PRINT_ASSERT(a,<=,3);
		return 10*a + Metric::index(i,j);
	}

	Christoffel(){
		data = NaN;
	}

	void contract2(const double kup[4], double result[4]) const{
		for(unsigned a=0; a<4; a++){
			result[a] = 0;
			for(unsigned i=0; i<4; i++)
				for(unsigned j=0; j<4; j++)
					result[a] += data[index(a,i,j)] * kup[i]*kup[j];
		}
	}
};
#endif
