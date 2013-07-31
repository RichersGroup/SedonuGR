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
}

//---------------------------------------------------------
// Initialize with passed vector
//---------------------------------------------------------
void locate_array::init(std::vector<double> a)
{
  x.resize(a.size());
  for (int i=0;i<x.size();i++) x[i] = a[i];
}

//---------------------------------------------------------
// locate
// If off the boundaries of the array, return the
// boundary value
//---------------------------------------------------------
int locate_array::locate(double z)
{
  // the degenerate case always returns 0
  if (x.size() == 1) return 0;
  
  // a form of locate from numerical recipes
  int bm;                             // mid point
  int bl = 0;                         // lower bound
  int bu = x.size()-1;               // upper bound

  // check if we are off the ends of the array
  if (z >= x[bu]) return (int)x.size();
  if (z <= x[bl]) return 0;
    
  // search the array for this index
  while (bu-bl > 1)
  {
    bm = (bu+bl)/2;
    if (x[bm] <= z) bl = bm;
    else bu = bm;
  }
  return bl;
} 


//---------------------------------------------------------
// Linear Interpolation of a passed array, find the zone
//---------------------------------------------------------
double locate_array::value_at(double z, std::vector<double> y)
{
  int ind = locate(z);
  return value_at(z, ind,y);
}

double locate_array::value_at(double z, int ind, std::vector<double> y)
{
  double v,slope;
  if (ind < x.size()-1)
  {
    int i1    = ind;
    int i2    = ind + 1;
    slope = (y[i2]-y[i1])/(x[i2] - x[i1]);
    v     = y[ind] + slope*(z - x[i1]);
  }
  else
  {
    int i2    = ind;
    int i1    = ind - 1;
    slope = (y[i2]-y[i1])/(x[i2] - x[i1]);
    v     = y[ind] + slope*(z - x[i1]);
  }

  return v;
}



//---------------------------------------------------------
// sample uniformally in zone
//---------------------------------------------------------
double locate_array::sample(int i, double z)
{
  if (i == x.size()) return x[i];
  return x[i] + (x[i+1] - x[i])*z;
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
  
