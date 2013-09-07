#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "locate_array.h"

using namespace std;


//---------------------------------------------------------
// Just allocation the memory for this
//---------------------------------------------------------
void locate_array::init(int n) 
{
  x.resize(n);
  for (int i=0;i<x.size();i++)  x[i] = 0;
}

//---------------------------------------------------------
// Initialize with start, stop and delta
//---------------------------------------------------------
void locate_array::init(double start, double stop, double del)
{
  int n = (stop - start)/del;
  if (n < 1) n = 1;
  x.resize(n);

  for (int i=0;i<x.size();i++) x[i] = start + i*del;
  do_log_interpolate = 0;
}

//---------------------------------------------------------
// Initialize with start, stop and n_pts
//---------------------------------------------------------
void locate_array::init(double start, double stop, int n)
{
  if (n < 1) n = 1;
  double del = (stop - start)/(1.0*n-1);  
  x.resize(n);

  for (int i=0;i<x.size();i++) x[i] = start + i*del;
  do_log_interpolate = 0;
}

//---------------------------------------------------------
// Initialize with passed vector
//---------------------------------------------------------
void locate_array::init(std::vector<double> a)
{
  x.resize(a.size());
  for (int i=0;i<x.size();i++) x[i] = a[i];
  do_log_interpolate = 0;
}

//---------------------------------------------------------
// locate (return closest index below the value)
// if off boundary, return -1 or size of array
//---------------------------------------------------------
int locate_array::locate(double xval)
{
  // a form of locate from numerical recipes
  int bm;                             // mid point
  int bl = 0;                         // lower bound
  int bu = x.size()-1;                // upper bound

  // check if we are off the ends of the array
  if (xval >= x[bu]) return x.size();
  if (xval <= x[bl]) return -1;
    
  // search the array for this index
  while (bu-bl > 1)
  {
    bm = (bu+bl)/2;
    if (x[bm] <= xval) bl = bm;
    else bu = bm;
  }
  return bl;
} 


//---------------------------------------------------------
// Linear Interpolation of a passed array, find the zone
//---------------------------------------------------------
double locate_array::interpolate_between(double xval, int i1, int i2, vector<double>& y)
{
  double slope = (y[i2]-y[i1]) / (x[i2]-x[i1]);
  double yval = y[i1] + slope*(xval - x[i1]);
  return yval;
}


//---------------------------------------------------------
// Log-Log Interpolation of a passed array, find the zone
//---------------------------------------------------------
double locate_array::log_interpolate_between(double xval, int i1, int i2, vector<double>& y)
{
  // safeguard against all opacities being 0
  // (e.g. if density is below nulib minimum)
  if(y[i1]==0 && y[i2]==0) 
    return 0;

  double slope = log(y[i2]/y[i1]) / log(x[i2]/x[i1]);
  double logyval = log(y[i1]) + slope*log(xval/x[i1]);

  return exp(logyval);
}


//---------------------------------------------------------
// sample uniformally in zone
//---------------------------------------------------------
double locate_array::sample(int i, double xval)
{
  if (i == x.size()) return x[i];
  return x[i] + (x[i+1] - x[i])*xval;
}

//---------------------------------------------------------
// simple printout
//---------------------------------------------------------
void locate_array::print()
{
  printf("# Print Locate Array; n_elements = %lu\n",x.size());
  for (int i=0;i<x.size();i++)
    printf("%4d %12.4e\n",i,x[i]);
}
  

 
//---------------------------------------------------------
// find the value of y at the locate_array's value of xval
// assumes 1-1 correspondence between y and locate_array
//---------------------------------------------------------
double locate_array::value_at(double xval, vector<double>& y){
  int ind = locate(xval);
  int i1, i2;
  if(ind < 0){                // If off left side of grid
    i1 = 0;
    i2 = 1;
  }
  else if(ind < x.size()-1){  // If within expected region of grid
    i1 = ind;
    i2 = ind + 1;
  }
  else{                       // If off the right side of the grid
    i1 = x.size() - 2;
    i2 = x.size() - 1;
  }

  if(do_log_interpolate) return log_interpolate_between(xval, i1, i2, y);
  else                   return     interpolate_between(xval, i1, i2, y);
}

