#ifndef _LOCATE_ARRAY_H
#define _LOCATE_ARRAY_H 

#include <vector>
#include "global_options.h"

class locate_array {

public:

	// where applicable, these values are the right bin wall (i.e. not the left)
	std::vector<double> x;
	double min;

	// other parameters
	InterpolationMethod interpolation_method;

	// constructors
	locate_array(int n=0);
	//locate_array(int n) {init(n);}

	// Return size of array (also, # of bins)
	unsigned size() const {return (int)x.size();}

	void init(const int size);
	void init(const double start,const double stop,const double del, const InterpolationMethod imeth=constant);
	void init(const double start,const double stop,const int n, const InterpolationMethod imeth=constant);
	void init(const std::vector<double> a, const double minval, const InterpolationMethod imeth=constant);
	void swap(locate_array new_array);

	// bottom of the bin left of i
	double bottom(const int i) const{
		return (i==0 ? min : x[i-1]);
	}

	// center of the bin left of i
	double center(const int i) const{
		return 0.5 * (bottom(i) + x[i]);
	}

	// width of the bin left of i
	double delta(const int i) const{
		return x[i] - bottom(i);
	}


	int    locate(const double) const;
	void   print() const;
	double value_at(const double nu, const std::vector<double>& array) const;

	// operators for easy access
	double  operator[] (const int i) const {return x[i];};
	double& operator[] (const int i)       {return x[i];};
	void resize(int i) {x.resize(i);};
};

#endif
