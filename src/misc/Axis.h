#ifndef _AXIS_H
#define _AXIS_H 1

#include <vector>
#include <algorithm>
#include "global_options.h"

using namespace std;

//======//
// AXIS //
//======//
class Axis{
public:
	double min;
	vector<double> top;
	vector<double> mid;

	Axis(const double min, vector<double>& top, vector<double>& mid){
		assert(top.size() == mid.size());
		this->min = min;
		this->top = top;
		this->mid = mid;
	}

	Axis() {
		min = NaN;
	}

	int size() const {return top.size();}

	int bin(const double x) const{
		if(x<min) return -1;
		else{
			// upper_bound returns first element greater than xval
			// values mark bin tops, so this is what we want
			int ind = upper_bound(top.begin(), top.end(), x) - top.begin();
			assert(ind>=0);
			assert(ind<=(int)top.size());
			return ind;
		}
	}
};

#endif
