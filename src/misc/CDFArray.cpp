/*
//  Copyright (c) 2015, California Institute of Technology and the Regents
//  of the University of California, based on research sponsored by the
//  United States Department of Energy. All rights reserved.
//
//  This file is part of Sedonu.
//
//  Sedonu is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  Neither the name of the California Institute of Technology (Caltech)
//  nor the University of California nor the names of its contributors 
//  may be used to endorse or promote products derived from this software
//  without specific prior written permission.
//
//  Sedonu is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with Sedonu.  If not, see <http://www.gnu.org/licenses/>.
//
*/

#include <cmath>
#include <algorithm>
#include <cstdio>
#include "global_options.h"
#include "CDFArray.h"

using namespace std;

// safe constructor
CDFArray::CDFArray(const int iorder){
	N = NaN;
	PRINT_ASSERT(iorder==1,||,iorder==3);
	interpolation_order = iorder;
}

//------------------------------------------------------
// return the actual y value, not the integrated
//------------------------------------------------------
double CDFArray::get_value(const int i) const
{
	PRINT_ASSERT(i,>=,0);
	PRINT_ASSERT(i,<,(int)size());
	if (i==0) return y[0];
	else return (y[i] - y[i-1]);
}





//------------------------------------------------------
// set the actual y value, not the integrated
// must be called in order
//------------------------------------------------------
void CDFArray::set_value(const int i, const double f)
{
	PRINT_ASSERT(i,>=,0);
	PRINT_ASSERT(i,<,(int)size());
	y[i] = ( i==0 ? f : y[i-1]+f );
}

//------------------------------------------------------
// Normalize such that the last entry is 1.0
//------------------------------------------------------
void CDFArray::normalize(double cutoff)
{
	// check for zero array, set to all constant
	if (y.back() == 0){
		y.assign(y.size(),1.0);
		N = 0;
	}
	// normalize to end = 1.0
	else{
		N = y.back();
		PRINT_ASSERT(N,>,0);
		double N_inv = 1.0/N;
		for(unsigned i=0;i<y.size();i++)   y[i] *= N_inv;
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
int CDFArray::get_index(const double yval) const
{
	PRINT_ASSERT(yval,>=,0);
	PRINT_ASSERT(yval,<=,1.0);
	PRINT_ASSERT(fabs(y.back()-1.0),<,1.0e-15);
	int i = upper_bound(y.begin(), y.end(), yval) - y.begin();
	PRINT_ASSERT(i,>=,0);
	PRINT_ASSERT(i,<=,(int)size());
	return i;
}

//-----------------------------------------------------
// return centered tangent at the corresponding vertex.
// i==-1 means at the left boundary
// assumes the cdf is monotonic
//-----------------------------------------------------
double CDFArray::inverse_tangent(const int i, const LocateArray* xgrid) const{
	int N = size();
	PRINT_ASSERT(i,>=,-1);
	PRINT_ASSERT(i,<=,(int)size()-1);

	// two-point stencil on boundary, 3-point stencil elsewhere. Keep out infinities.
	double result = 0.0;
	double left_secant  = (i==-1  ? numeric_limits<double>::infinity() : inverse_secant(i-1,i,xgrid));
	double right_secant = (i==N-1 ? numeric_limits<double>::infinity() : inverse_secant(i,i+1,xgrid));
	PRINT_ASSERT(!isinf<bool>(left_secant),||,!isinf<bool>(right_secant));
	if     (i==-1 ) result = right_secant;
	else if(i==N-1) result = left_secant;
	else{
		if     (isinf<bool>(left_secant )) result = right_secant;
		else if(isinf<bool>(right_secant)) result = left_secant;
		else                         result = 0.5 * (right_secant + left_secant);
	}
	return result;
}
double CDFArray::tangent(const int i, const LocateArray* xgrid) const{
	int N = size();
	PRINT_ASSERT(i,>=,-1);
	PRINT_ASSERT(i,<=,(int)size()-1);

	// two-point stencil on boundary, 3-point stencil elsewhere. Keep out infinities.
	double result = 0.0;
	double left_secant  = (i==-1  ? numeric_limits<double>::infinity() : secant(i-1,i,xgrid));
	double right_secant = (i==N-1 ? numeric_limits<double>::infinity() : secant(i,i+1,xgrid));
	PRINT_ASSERT(!isinf<bool>(left_secant),||,!isinf<bool>(right_secant));
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
double CDFArray::inverse_secant(const int i, const int j, const LocateArray* xgrid) const{
	return 1.0 / secant(i,j,xgrid);
}
double CDFArray::secant(const int i, const int j, const LocateArray* xgrid) const{
	PRINT_ASSERT(i>=-1,&&,i<(int)size()-1);
	PRINT_ASSERT(j>=0, &&,j<(int)size()  );
	PRINT_ASSERT(j,>,i);

	double result = 0.0;
	if(i==-1) result = (y[j]-0.0 ) / ((*xgrid)[j]-xgrid->min );
	else      result = (y[j]-y[i]) / ((*xgrid)[j]-(*xgrid)[i]);
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
	PRINT_ASSERT(0.0,<=,t);
	PRINT_ASSERT(t,<=,1.0);
	return (1.0+2.0*t)*(1.0-t)*(1.0-t);
}
double h10(const double t){
	PRINT_ASSERT(0,<=,t);
	PRINT_ASSERT(t,<=,1.0);
	return t*(1.0-t)*(1.0-t);
}
double h01(const double t){
	PRINT_ASSERT(0.0,<=,t);
	PRINT_ASSERT(t,<=,1.0);
	return t*t*(3.0-2.0*t);
}
double h11(const double t){
	PRINT_ASSERT(0.0,<=,t);
	PRINT_ASSERT(t,<=,1.0);
	return t*t*(t-1.0);
}
double h00p(const double t){
	PRINT_ASSERT(0.0,<=,t);
	PRINT_ASSERT(t,<=,1.0);
	return 6.0 * t * (t-1.0);
}
double h10p(const double t){
	PRINT_ASSERT(0,<=,t);
	PRINT_ASSERT(t,<=,1.0);
	return (1.0-t)*(3.0-t);
}
double h01p(const double t){
	PRINT_ASSERT(0.0,<=,t);
	PRINT_ASSERT(t,<=,1.0);
	return 6.0*t*(t-1.0);
}
double h11p(const double t){
	PRINT_ASSERT(0.0,<=,t);
	PRINT_ASSERT(t,<=,1.0);
	return t*(3.0*t-2.0);
}
double CDFArray::invert(const double rand, const LocateArray* xgrid, const int i_in) const{
	double result = 0;
	assert(interpolation_order==1 || interpolation_order==3 || interpolation_order==0);
	if     (interpolation_order==0) result = invert_piecewise(rand,xgrid,i_in);
	else if(interpolation_order==1) result = invert_linear(rand,xgrid,i_in);
	else if(interpolation_order==3) result = invert_cubic(rand,xgrid,i_in);
	PRINT_ASSERT(result,>,0);
	return result;
}
double CDFArray::interpolate_pdf(const double x, const LocateArray* xgrid) const
{
	double result = 0;
	assert(interpolation_order==1 || interpolation_order==3 || interpolation_order==0);
	if     (interpolation_order==0) result = interpolate_pdf_piecewise(x,xgrid);
	else if(interpolation_order==1) result = interpolate_pdf_linear(x,xgrid);
	else if(interpolation_order==3) result = interpolate_pdf_cubic(x,xgrid);
	PRINT_ASSERT(result,>=,0);
	return result;
}
double CDFArray::invert_cubic(const double rand, const LocateArray* xgrid, const int i_in) const
// INCONSISTENCY - the emissivity is integrated assuming piecewise constant.
// This is inconsistent with cubic interpolation.
{
	PRINT_ASSERT(rand,>=,0);
	PRINT_ASSERT(rand,<=,1);
	int i = (i_in<0 ? get_index(rand) : i_in);
	PRINT_ASSERT(i,<,(int)size());
	PRINT_ASSERT(i,>=,0);

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
	double slope_inv = secant(i-1,i,xgrid);
	double slope = 1.0/slope_inv;
	double limiter = 3.0;
	double alpha = mLeft*slope_inv;
	double beta = mRight*slope_inv;
	PRINT_ASSERT(alpha,>,0);
	PRINT_ASSERT(beta,>,0);
	if(alpha*alpha + beta*beta > limiter*limiter){
		double tau = limiter/sqrt(alpha*alpha+beta*beta);
		mLeft = tau*alpha*slope;
		mRight = tau*beta*slope;
	}

	// return interpolated function
	double h = yRight-yLeft;
	double h_inv = 1.0/h;
	double t = 0;
	if(i_in<0) t = (rand-yLeft)*h_inv;
	else t = rand*(yRight-yLeft)*h_inv;
	PRINT_ASSERT(t,>=,0);
	PRINT_ASSERT(t,<=,1);
	double result = xLeft*h00(t) + h*mLeft*h10(t) + xRight*h01(t) + h*mRight*h11(t);
	result = max(xLeft,result);
	result = min(xRight,result);
	PRINT_ASSERT(xLeft,<=,result);
	PRINT_ASSERT(xRight,>=,result);
	return result;
}
double CDFArray::interpolate_pdf_cubic(const double x, const LocateArray* xgrid) const
{
	PRINT_ASSERT(x,>=,xgrid->min);
	PRINT_ASSERT(x,<=,xgrid->x[xgrid->size()-1]);

	// get the upper index
	int i = xgrid->locate(x);
	PRINT_ASSERT(i,<,(int)size());
	PRINT_ASSERT(i,>=,0);

	// check for degenerate case (left and right values are equal)
	double yRight = y[i];
	double xRight = (*xgrid)[i];
	double yLeft = (i>0 ?        y[i-1] : 0         );
	double xLeft = (i>0 ? (*xgrid)[i-1] : xgrid->min);
	if(yRight == yLeft) return yRight;

	// get left and right tangents
	double mLeft  = tangent(i-1,xgrid);
	double mRight = tangent(i  ,xgrid);
	assert(!isinf<bool>(mLeft ));
	assert(!isinf<bool>(mRight));

	// prevent overshoot, ensure monotonicity
	double slope_inv = secant(i-1,i,xgrid);
	double slope = 1.0/slope_inv;
	double limiter = 3.0;
	double alpha = mLeft*slope_inv;
	double beta = mRight*slope_inv;
	PRINT_ASSERT(alpha,>,0);
	PRINT_ASSERT(beta,>,0);
	if(alpha*alpha + beta*beta > limiter*limiter){
		double tau = limiter/sqrt(alpha*alpha+beta*beta);
		mLeft = tau*alpha*slope;
		mRight = tau*beta*slope;
	}

	// return interpolated function
	double h = xRight-xLeft;
	double t = (x-xLeft)/h;
	double dtdx = h;
	PRINT_ASSERT(t,>=,0);
	PRINT_ASSERT(t,<=,1);
	double result = yLeft*h00(t) + h*mLeft*h10(t) + yRight*h01(t) + h*mRight*h11(t);
	//double result = dtdx * (yLeft*h00p(t) + h*mLeft*h10p(t) + yRight*h01p(t) + h*mRight*h11p(t));
	result = max(yLeft,result);
	result = min(yRight,result);
	//PRINT_ASSERT(yLeft,<=,result);
	//PRINT_ASSERT(yRight,>=,result);
	return result;
}
double CDFArray::interpolate_pdf_linear(const double x, const LocateArray* xgrid) const
{
	PRINT_ASSERT(x,>=,xgrid->min);
	PRINT_ASSERT(x,<=,xgrid->x[xgrid->size()-1]);

	// get the upper/lower indices
	int upper = xgrid->locate(x);
	PRINT_ASSERT(upper,>=,0);
	int lower = upper-1;
	PRINT_ASSERT(lower,<,(int)xgrid->size());

	// get the x and y values of the left and right sides
	double x1,x2,y1,y2;
	if(upper==0){
		x1 = xgrid->min;
		x2 = (*xgrid)[0];
		y1 = 0;
		y2 = y[0];
	}
	else{
		x1 = (*xgrid)[lower];
		x2 = (*xgrid)[upper];
		y1 = y[lower];
		y2 = y[upper];
	}

	double slope = (y2-y1)/(x2-x1);
	double result = slope; //y1 + slope*(x-x1);
	//PRINT_ASSERT(result <= (*xgrid)[xgrid->size()-1]);
	//PRINT_ASSERT(result >= xgrid->min);
	return result;
}
double CDFArray::invert_linear(const double rand, const LocateArray* xgrid, const int i_in) const
{
	PRINT_ASSERT(rand,>=,0);
	PRINT_ASSERT(rand,<=,1);
	int i = (i_in<0 ? get_index(rand) : i_in);
	PRINT_ASSERT(i,<,(int)size());
	PRINT_ASSERT(i,>=,0);

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
	double result = 0;
	if(i_in<0) result = xLeft + slope*(rand-yLeft);
	else result = xLeft + rand*slope*(yRight-yLeft);
	PRINT_ASSERT(xLeft,<=,result);
	PRINT_ASSERT(xRight,>=,result);
	return result;
}
double CDFArray::invert_piecewise(const double rand, const LocateArray* xgrid, const int i_in) const
{
	PRINT_ASSERT(rand,>=,0);
	PRINT_ASSERT(rand,<=,1);
	int i = (i_in<0 ? get_index(rand) : i_in);
	PRINT_ASSERT(i,<,(int)size());
	PRINT_ASSERT(i,>=,0);

	return xgrid->center(i);
}
double CDFArray::interpolate_pdf_piecewise(const double x, const LocateArray* xgrid) const
{
	PRINT_ASSERT(x,>=,xgrid->min);
	PRINT_ASSERT(x,<=,xgrid->x[xgrid->size()-1]);

	//if(x<(*xgrid)[0]) return 0.0;

	int upper = xgrid->locate(x);
	PRINT_ASSERT(upper,>=,0);
	int lower = upper-1;
	PRINT_ASSERT(lower,<,(int)xgrid->size());

	return get_value(upper); //y[lower];
}

//------------------------------------------------------
// Simple printout
//------------------------------------------------------
void CDFArray::print() const{
	for(unsigned i=0;i<y.size();i++)
		printf("%5d %10.4e %10.4e\n",i,get_value(i),y[i]);
}

//------------------------------------------------------
// Clear the arrays
//------------------------------------------------------
void CDFArray::wipe()
{
	y.assign(y.size(), 1.0);
}

//------------------------------------------------------------
// just returning the size of the array
//------------------------------------------------------------
unsigned CDFArray::size() const
{
	return y.size();
}
