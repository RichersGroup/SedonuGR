#ifndef _LOCATE_ARRAY_H
#define _LOCATE_ARRAY_H 

#include <vector>
#include "cdf_array.h"


class locate_array {

public:

  std::vector<double> x;
  int do_log_interpolate;

  // constructors
  locate_array()  {}
  locate_array(int n) {init(n);}

  // Return
  int    size()        {return (int)x.size();}

  void init(int);
  void init(double,double,double);
  void init(double,double,int);
  void init(std::vector<double>);

  double center(int i)   {
    if (i >= x.size()-1) {return x.back();}
    else return 0.5*(x[i] + x[i+1]);}

  double delta(int i)   {
    if (i >= x.size()-1) {return x[i] - x[i-1]; }
    else return (x[i+1] - x[i]); }


  int    locate(double);
  double interpolate_between(double,int,int,std::vector<double>&);
  double log_interpolate_between(double,int,int,std::vector<double>&);
  double sample(int, double);
  void   print();
  double value_at(double nu, std::vector<double>& array);

  // operators for easy access
  double& operator[] (const int i) {return x[i];};
  void resize(int i) {x.resize(i);};
};

#endif
