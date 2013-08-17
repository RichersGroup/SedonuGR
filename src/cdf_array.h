#ifndef _CDF_H
#define _CDF_H 1

#include <vector>

//**********************************************************
// CDF == Comulative Distribution Function
//
// This simple class just holds a vector which should be
// monitonically increasing and reaches unity
// We can sample from it using a binary search.
//**********************************************************

class cdf_array
{

private:
  
  std::vector<double> y;
  double N;
  
public:

  void resize(int n)  {y.resize(n); }

  double get(int i)            {return y[i];}   // Get local CDF value
  void   set(int i, double f)  {y[i] = f;}      // Set cell CDF value 

  void   set_value(int i, double f);     // set the actual (not CDF) value
  double get_value(int i);                   // Get the actual (not CDF) value
 
  void normalize();         // normalize the cdf, so that final value = 1
  double get_N();           // get the normalization factor
  int  sample(double z);    // sample from the CDF, when passed a random #
  void print();             
  void wipe();
  void MPI_combine(); 
  int size();

};

#endif
