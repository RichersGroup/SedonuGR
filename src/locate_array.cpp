#include "global_options.h"
#include <vector>
#include <algorithm>
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "locate_array.h"


//-----------------------------
// safe nan-filled constructor
//-----------------------------
locate_array::locate_array(const int n){
	min = NaN;
	interpolation_method = constant;
	x.assign(n,NaN);
}

//---------------------------------------------------------
// Initialize with start, stop and delta
// if start==stop make it a catch-all
//---------------------------------------------------------
void locate_array::init(const double start, const double stop, const double del, const InterpolationMethod imeth)
{
	assert(stop >= start);
	assert(del > 0);
	if(start==stop){
		x.resize(1);
		min = -numeric_limits<double>::infinity();
		x[0] = numeric_limits<double>::infinity();
	}

	else{
		int n = ceil( (stop-start)/del );
		n = max(n,1);
		x.resize(n);
		interpolation_method = imeth;

		min = start;
        #pragma omp parallel for
		for(int i=0; i<n-1; i++){
			x[i] = start + (double)(i+1)*del;
			if(i>0) assert(x[i] > x[i-1]);
		}
		x[n-1] = stop;
		assert(x[n-1] > x[n-2]);
	}
}

//---------------------------------------------------------
// Initialize with start, stop and n_pts
// if start==stop make it a catch-all
// if n==0 make it a catch-all
//---------------------------------------------------------
void locate_array::init(const double start, const double stop, const int n, const InterpolationMethod imeth)
{
	assert(stop >= start);
	assert(n >= 0);
	if(start==stop || n==0){
		x.resize(1);
		min = -numeric_limits<double>::infinity();
		x[0] = numeric_limits<double>::infinity();
	}

	else{
		double del = (stop - start)/(double)n;
		x.resize(n);
		interpolation_method = constant;

		min = start;
        #pragma omp parallel for
		for(int i=0; i<n-1; i++) x[i] = start + (i+1)*del;
		x[n-1] = stop;
		assert(x[n-1] > x[n-2]);
	}
}

//---------------------------------------------------------
// Initialize with passed vector
//---------------------------------------------------------
void locate_array::init(const std::vector<double> a, const double minval, const InterpolationMethod imeth)
{
	min = minval;
	interpolation_method = constant;
	x.assign(a.begin(), a.end());
}

//---------------------------------------------------------
// locate (return closest index below the value)
// if off left side of boundary, returns 0
// if off right side of boundary, returns size
//---------------------------------------------------------
int locate_array::locate(const double xval) const
{
	// upper_bound returns first element greater than xval
	// values mark bin tops, so this is what we want
	int ind = upper_bound(x.begin(), x.end(), xval) - x.begin();
	assert(ind >= 0);
	assert(ind <= (int)x.size());
	return ind;
} 


//---------------------------------------------------------
// simple printout
//---------------------------------------------------------
void locate_array::print() const
{
	printf("# Print Locate Array; n_elements = %lu\n",x.size());
	printf("min %12.4e\n",min);
	for(unsigned i=0;i<x.size();i++)
		printf("%4d %12.4e\n",i,x[i]);
}



//---------------------------------------------------------
// find the value of y at the locate_array's value of xval
// assumes y array values are stored at center of locate_array bins
//---------------------------------------------------------
double lin_interpolate_between(const double xval, const double xleft, const double xright, const double yleft, const double yright){
	assert(xleft<xright);
	double slope = (yright-yleft) / (xright-xleft);
	double yval = yleft + slope*(xval - xleft);
	return yval;
}
double log_interpolate_between(const double xval, const double xleft, const double xright, const double yleft, const double yright){
	assert(xleft<xright);
	assert(xleft  > 0);
	assert(xright > 0);
	assert(yleft  > 0);
	assert(yright > 0);

	// safeguard against equal opacities
	if(yleft==yright) return yleft;

	// do logarithmic interpolation
	double slope = log(yright/yleft) / log(xright/xleft);
	double logyval = log(yleft) + slope*log(xval/xleft);
	return exp(logyval);
}
double pow_interpolate_between(const double xval, const double xleft, const double xright, const double yleft, const double yright){
	assert(xleft<xright);
	assert(xleft  > 0);
	assert(xright > 0);

	// safeguard against equal opacities
	if(yleft==yright) return yleft;

	// do power interpolation
	double exponent = (xval-xleft)/(xright-xleft);
	double base = (yright/yleft);
	return yleft * pow(base,exponent);
}
double locate_array::value_at(const double xval, const vector<double>& y) const{
	int ind = locate(xval);
	assert(ind >= 0);
	assert(ind <= (int)x.size());

	if(interpolation_method == constant) return y[ind];
	else{
		int i1, i2;
		if(ind == 0 || xval<=center(0)){   // If off left side of grid
			i1 = 0;
			i2 = 1;
		}
		else if(ind==(int)x.size() || xval>=center(x.size()-1)){ // If off the right side of the grid
			i1 = x.size() - 2;
			i2 = x.size() - 1;
		}
		else{    // If within expected region of grid
			if(xval<center(ind)){
				i1 = ind - 1;
				i2 = ind;
			}
			else{
				i1 = ind;
				i2 = ind + 1;
			}
		}
		double xleft  = center(i1);
		double xright = center(i2);
		double yleft  = y[i1];
		double yright = y[i2];


		if     (interpolation_method == linear     ) return lin_interpolate_between(xval,xleft,xright,yleft,yright);
		else if(interpolation_method == logarithmic) return log_interpolate_between(xval,xleft,xright,yleft,yright);
		else{ assert(interpolation_method == power );return pow_interpolate_between(xval,xleft,xright,yleft,yright);}
	}
}


void locate_array::swap(locate_array& new_array){
	// swap the vectors
	x.swap(new_array.x);

	// swap the minimum values
	double min_tmp = min;
	min = new_array.min;
	new_array.min = min_tmp;

	// swap the do_log_interpolate parameters
	InterpolationMethod tmp = interpolation_method;
	interpolation_method = new_array.interpolation_method;
	new_array.interpolation_method = tmp;
}

void locate_array::copy(const locate_array& new_array){
	// copy the vector
	x.resize(new_array.x.size());
	for(unsigned i=0; i<new_array.x.size(); i++) x[i] = new_array.x[i];

	// copy the minimum value
	min = new_array.min;

	// copy the do_log_interpolate parameter
	interpolation_method = new_array.interpolation_method;
}
