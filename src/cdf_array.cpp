#include <mpi.h>
#include "cdf_array.h"


//------------------------------------------------------
// return the actual y value, not the integrated
//------------------------------------------------------
double cdf_array::get_value(int i)   
{
  if (i==0) return y[0];
  else return (y[i] - y[i-1]);  
}



//------------------------------------------------------
// set the actual y value, not the integrated
//------------------------------------------------------
void cdf_array::set_value(int i, double f)   
{
  if (i==0) y[0] = f;
  else y[i] = y[i-1] + f;
}


//------------------------------------------------------
// Normalize such that the last entry is 1.0
//------------------------------------------------------
void cdf_array::normalize() 
{
  // check for zero array, set to all constant
  if (y.back() == 0) 
    for (int i=1;i<y.size();i++) y[i] = y[i-1] + 1;

  // normalize to end = 1.0
  N = y.back();
  for (int i=0;i<y.size();i++)   y[i] /= N;
}


//---------------------------------------------------------
// Sample the probability distribution using binary search.
// Pass a random number betwen 0 and 1.  
// Returns the bin number
//---------------------------------------------------------
int cdf_array::sample(double z)
{
  // mid, lower, and upper points
  int bm;                         // mid point
  int bl = 0;                     // lower bound
  int bu = y.size() - 1;        // upper bound
  // see if we are off the top or bottom
  if (z > y[bu]) return bu;
  if (z < y[bl]) return bl;
  // search
  while (bu-bl > 1)
  {
    bm = (bu+bl)/2;
    if (y[bm] <= z) bl = bm;
    else bu = bm;
  }
  return bu;
}


//------------------------------------------------------
// Simple printout
//------------------------------------------------------
void cdf_array::print() {
  for (int i=0;i<y.size();i++) 
    printf("%5d %10.4e %10.4e\n",i,get_value(i),y[i]);
}
  
//------------------------------------------------------
// Clear the arrays
//------------------------------------------------------
void cdf_array::wipe()
{
  for (int i=0;i<y.size();i++)  y[i] = 0;
}
  
//------------------------------------------------------
// MPI Reduce this class
//------------------------------------------------------
void cdf_array::MPI_combine() 
{
  // double *new_ptr;
  // int i;          
    
  // // allocate the memory for new pointer
  // new_ptr = new Type[n_elements];
  // // zero out array
  // for (i=0;i<n_elements;i++) new_ptr[i] = 0;
  // MPI_Allreduce(array,new_ptr,n_elements,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  // // put back into place
  // for (i=0;i<n_elements;i++) array[i] = new_ptr[i];
  // // free up the memory
  // delete new_ptr;
}   


//------------------------------------------------------------
// return the normalization variable
// should not be set by anything but normalize()
//------------------------------------------------------------
double cdf_array::get_N()
{
  return N;
}
