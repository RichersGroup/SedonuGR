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

#include <algorithm>
#include <cmath>
#include <cstdio>
#include "LocateArray.h"
#include "global_options.h"

using namespace std;

//-----------------------------
// safe nan-filled constructor
//-----------------------------
LocateArray::LocateArray(const int n){
	min = NaN;
	interpolation_method = constant;
	x.assign(n,NaN);
}

//---------------------------------------------------------
// Initialize with start, stop and delta
// if start==stop make it a catch-all
//---------------------------------------------------------
void LocateArray::init(const double start, const double stop, const double del, const InterpolationMethod imeth)
{
	PRINT_ASSERT(stop,>=,start);
	PRINT_ASSERT(del,>,0);
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
			if(i>0) PRINT_ASSERT(x[i],>,x[i-1]);
		}
		x[n-1] = stop;
		PRINT_ASSERT(x[n-1],>,x[n-2]);
	}
}

//---------------------------------------------------------
// Initialize with start, stop and n_pts
// if start==stop make it a catch-all
// if n==0 make it a catch-all
//---------------------------------------------------------
void LocateArray::init(const double start, const double stop, const int n, const InterpolationMethod imeth)
{
	PRINT_ASSERT(stop,>=,start);
	PRINT_ASSERT(n,>=,0);
	if(start==stop || n==0){
		x.resize(1);
		min = -numeric_limits<double>::infinity();
		x[0] = numeric_limits<double>::infinity();
	}

	else{
		double del = (stop - start)/(double)n;
		x.resize(n);
		interpolation_method = imeth;

		min = start;
        #pragma omp parallel for
		for(int i=0; i<n-1; i++) x[i] = start + (i+1)*del;
		x[n-1] = stop;
		PRINT_ASSERT(x[n-1],>,x[n-2]);
	}
}

//---------------------------------------------------------
// Initialize with passed vector
//---------------------------------------------------------
void LocateArray::init(const std::vector<double> a, const double minval, const InterpolationMethod imeth)
{
	min = minval;
	interpolation_method = imeth;
	x.assign(a.begin(), a.end());
}

//---------------------------------------------------------
// locate (return closest index below the value)
// if off left side of boundary, returns 0
// if off right side of boundary, returns size
//---------------------------------------------------------
int LocateArray::locate(const double xval) const
{
	// upper_bound returns first element greater than xval
	// values mark bin tops, so this is what we want
	int ind = upper_bound(x.begin(), x.end(), xval) - x.begin();
	PRINT_ASSERT(ind,>=,0);
	PRINT_ASSERT(ind,<=,(int)x.size());
	return ind;
} 


//---------------------------------------------------------
// simple printout
//---------------------------------------------------------
void LocateArray::print() const
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
	PRINT_ASSERT(xleft,<,xright);
	double slope = (yright-yleft) / (xright-xleft);
	double yval = yleft + slope*(xval - xleft);
	return yval;
}
double log_interpolate_between(const double xval, const double xleft, const double xright, const double yleft, const double yright){
	PRINT_ASSERT(xleft,<,xright);
	PRINT_ASSERT(xleft,>,0);
	PRINT_ASSERT(xright,>,0);
	PRINT_ASSERT(yleft,>,0);
	PRINT_ASSERT(yright,>,0);

	// safeguard against equal opacities
	if(yleft==yright) return yleft;

	// do logarithmic interpolation
	double slope = log(yright/yleft) / log(xright/xleft);
	double logyval = log(yleft) + slope*log(xval/xleft);
	return exp(logyval);
}
double pow_interpolate_between(const double xval, const double xleft, const double xright, const double yleft, const double yright){
	PRINT_ASSERT(xleft,<,xright);
	PRINT_ASSERT(xleft,>,0);
	PRINT_ASSERT(xright,>,0);

	// safeguard against equal opacities
	if(yleft==yright) return yleft;

	// do power interpolation
	double exponent = (xval-xleft)/(xright-xleft);
	double base = (yright/yleft);
	double retval = yleft * pow(base,exponent);
	return retval;
}
double LocateArray::value_at(const double xval, const vector<double>& y) const{
	int ind = locate(xval);
	PRINT_ASSERT(ind,>=,0);
	PRINT_ASSERT(ind,<=,(int)x.size());

	double retval = -1;
	if(interpolation_method == constant){
		if(ind==x.size()) ind -= 1;
		retval = y[ind];
	}
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

		if     (interpolation_method == linear     ) retval = lin_interpolate_between(xval,xleft,xright,yleft,yright);
		else if(interpolation_method == logarithmic) retval = log_interpolate_between(xval,xleft,xright,yleft,yright);
		else{ PRINT_ASSERT(interpolation_method,==,power );retval = pow_interpolate_between(xval,xleft,xright,yleft,yright);}
	}
	PRINT_ASSERT(abs(retval),<,INFINITY);
	return retval;
}


void LocateArray::swap(LocateArray& new_array){
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

void LocateArray::copy(const LocateArray& new_array){
	// copy the vector
	x.resize(new_array.x.size());
	for(unsigned i=0; i<new_array.x.size(); i++) x[i] = new_array.x[i];

	// copy the minimum value
	min = new_array.min;

	// copy the do_log_interpolate parameter
	interpolation_method = new_array.interpolation_method;
}
