#ifndef _LOCATE_ARRAY_H
#define _LOCATE_ARRAY_H 

#include <vector>

class locate_array {

public:

  // where applicable, these values are the right bin wall (i.e. not the left)
  std::vector<double> x;
  double min;

  // other parameters
  int do_log_interpolate;

  // constructors
  locate_array()  {}
  locate_array(int n) {init(n);}

  // Return size of array (also, # of bins)
  int size() {return (int)x.size();}

  void init(int);
  void init(double,double,double);
  void init(double,double,int);
  void init(std::vector<double>, double minval);
  void swap(locate_array new_array);

  // center of the bin left of nu_i
  double center(int i){
    if (i == 0) return 0.5*(min    + x[0]);
    else        return 0.5*(x[i-1] + x[i]);
  }

  // width of the bin left of nu_i
  double delta(int i){
    if (i == 0) return x[0] - min;
    else        return x[i] - x[i-1];
  }


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
