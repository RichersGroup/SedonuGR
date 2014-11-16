#include "global_options.h"
#include <algorithm>
#include <stdio.h>
#include <cmath>
#include "cdf_array.h"
#include "locate_array.h"


// safe constructor
cdf_array::cdf_array(const int iorder){
	N = NaN;
	assert(iorder==1 || iorder==3);
	interpolation_order = iorder;
}

//------------------------------------------------------
// return the actual y value, not the integrated
//------------------------------------------------------
double cdf_array::get_value(const int i) const
{
	assert(i >= 0);
	assert(i < (int)size());
	if (i==0) return y[0];
	else return (y[i] - y[i-1]);
}

//------------------------------------------------------
// set the actual y value, not the integrated
// must be called in order
//------------------------------------------------------
void cdf_array::set_value(const int i, const double f)
{
	assert(i >= 0);
	assert(i < (int)size());
	y[i] = ( i==0 ? f : y[i-1]+f );
}

//------------------------------------------------------
// Normalize such that the last entry is 1.0
//------------------------------------------------------
void cdf_array::normalize(double cutoff)
{
	// check for zero array, set to all constant
	if (y.back() == 0){
		y.assign(y.size(),1.0);
		N = 0;
	}
	// normalize to end = 1.0
	else{
		N = y.back();
		assert(N > 0);
		for(unsigned i=0;i<y.size();i++)   y[i] /= ( N>0 ? N : 1 );
	}

	// set values below cutoff to 0
	if(cutoff>0){
		for(unsigned i=0; i<y.size(); i++){
			double val = get_value(i);
			if(val < cutoff) set_value(i,0);
			else set_value(i,val);
		}
		double tmp=N;
		normalize();
		N *= tmp;
	}
}

//---------------------------------------------------------
// Sample the probability distribution using binary search.
// Pass a number betwen 0 and 1.
// Returns the index of the first value larger than yval
// if larger than largest element, returns size
//---------------------------------------------------------
int cdf_array::get_index(const double yval) const
{
	assert(yval >= 0);
	assert(yval <= 1.0);
	int i = upper_bound(y.begin(), y.end(), yval) - y.begin();
	assert(i >= 0);
	assert(i < (int)size());
	return i;
}

//-----------------------------------------------------
// return centered tangent at the corresponding vertex.
// i==-1 means at the left boundary
// assumes the cdf is monotonic
//-----------------------------------------------------
double cdf_array::inverse_tangent(const int i, const locate_array* xgrid) const{
	int N = size();
	assert(i>=-1);
	assert(i<=(int)size()-1);

	// two-point stencil on boundary, 3-point stencil elsewhere. Keep out infinities.
	double result = 0.0;
	double left_secant  = (i==-1  ? numeric_limits<double>::infinity() : inverse_secant(i-1,i,xgrid));
	double right_secant = (i==N-1 ? numeric_limits<double>::infinity() : inverse_secant(i,i+1,xgrid));
	assert(!isinf<bool>(left_secant) || !isinf<bool>(right_secant));
	if     (i==-1 ) result = right_secant;
	else if(i==N-1) result = left_secant;
	else{
		if     (isinf<bool>(left_secant )) result = right_secant;
		else if(isinf<bool>(right_secant)) result = left_secant;
		else                         result = 0.5 * (right_secant + left_secant);
	}
	return result;
}

//--------------------------------------
// return secant line between two points
// (CDF is x value, xgrid is y value)
//--------------------------------------
double cdf_array::inverse_secant(const int i, const int j, const locate_array* xgrid) const{
	assert(i>=-1 && i<(int)size()-1);
	assert(j>=0  && j<(int)size()  );
	assert(j>i);

	double result = 0.0;
	if(i==-1) result = ((*xgrid)[j]-xgrid->min ) / (y[j]-0.0 );
	else      result = ((*xgrid)[j]-(*xgrid)[i]) / (y[j]-y[i]);
	return result;
}

//---------------------------------------------------------
// Sample the probability distribution using binary search.
// Pass a random number betwen 0 and 1.
// Returns the corresponding cubic-interpolated xgrid value
// assumes the CDF is monotonic
//---------------------------------------------------------
// first, the cubic hermit spline functions
double h00(const double t){
	assert(0.0<=t && t<=1.0);
	return (1.0+2.0*t)*(1.0-t)*(1.0-t);
}
double h10(const double t){
	assert(0<=t && t<=1.0);
	return t*(1.0-t)*(1.0-t);
}
double h01(const double t){
	assert(0.0<=t && t<=1.0);
	return t*t*(3.0-2.0*t);
}
double h11(const double t){
	assert(0.0<=t && t<=1.0);
	return t*t*(t-1.0);
}
double cdf_array::invert_cubic(const double rand, const locate_array* xgrid, const int i_in) const
// INCONSISTENCY - the emissivity is integrated assuming piecewise constant.
// This is inconsistent with cubic interpolation.
{
	assert(rand>=0 && rand<=1);
	int i = (i_in<0 ? get_index(rand) : i_in);
	assert(i<(int)size());
	assert(i>=0);

	// check for degenerate case (left and right values are equal)
	double yRight = y[i];
	double xRight = (*xgrid)[i];
	double yLeft = (i>0 ?        y[i-1] : 0         );
	double xLeft = (i>0 ? (*xgrid)[i-1] : xgrid->min);
	if(yRight == yLeft) return yRight;

	// get left and right tangents
	double mLeft  = inverse_tangent(i-1,xgrid);
	double mRight = inverse_tangent(i  ,xgrid);
	assert(!isinf<bool>(mLeft ));
	assert(!isinf<bool>(mRight));

	// prevent overshoot, ensure monotonicity
	double slope = inverse_secant(i-1,i,xgrid);
	double limiter = 3.0;
	double alpha = mLeft/slope;
	double beta = mRight/slope;
	assert(alpha>0);
	assert(beta>0);
	if(alpha*alpha + beta*beta > limiter*limiter){
		double tau = limiter/sqrt(alpha*alpha+beta*beta);
		mLeft = tau*alpha*slope;
		mRight = tau*beta*slope;
	}

	// return interpolated function
	double h = y[i]-yLeft;
	double t = (rand-yLeft)/h;
	assert(t>=0 && t<=1);
	double result = xLeft*h00(t) + h*mLeft*h10(t) + xRight*h01(t) + h*mRight*h11(t);
	result = max(xLeft,result);
	result = min(xRight,result);
	assert(xLeft<=result);
	assert(xRight>=result);
	return result;
}
double cdf_array::invert_linear(const double rand, const locate_array* xgrid, const int i_in) const
{
	assert(rand>=0 && rand<=1);
	int i = (i_in<0 ? get_index(rand) : i_in);
	assert(i<(int)size());
	assert(i>=0);

	// check for degenerate case (left and right values are equal)
	double xLeft = (i>0 ? (*xgrid)[i-1] : xgrid->min);
	double yLeft = (i>0 ?        y[i-1] : 0         );
	double yRight = y[i];
	double xRight = (*xgrid)[i];
	if(yRight == yLeft) return yRight;

	// get the slope between the two adjacent points
	double slope = inverse_secant(i-1,i,xgrid);
	assert(!isinf<bool>(slope));

	// return interpolated function
	double result = xLeft + slope*(rand-yLeft);
	assert(xLeft<=result);
	assert(xRight>=result);
	return result;
}

//------------------------------------------------------
// Simple printout
//------------------------------------------------------
void cdf_array::print() const{
	for(unsigned i=0;i<y.size();i++)
		printf("%5d %10.4e %10.4e\n",i,get_value(i),y[i]);
}

//------------------------------------------------------
// Clear the arrays
//------------------------------------------------------
void cdf_array::wipe()
{
	y.assign(y.size(), 1.0);
}

//------------------------------------------------------------
// just returning the size of the array
//------------------------------------------------------------
unsigned cdf_array::size() const
{
	return y.size();
}
