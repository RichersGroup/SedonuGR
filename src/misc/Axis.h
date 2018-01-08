#ifndef _AXIS_H
#define _AXIS_H 1

#include <vector>
#include <algorithm>
#include "global_options.h"
#include "hdf5.h"
#include <string>

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
		PRINT_ASSERT(top.size(),==,mid.size());
		this->min = min;
		this->top = top;
		this->mid = mid;
		for(int i=0; i<top.size(); i++){
			PRINT_ASSERT(top[i],>,mid[i]);
			PRINT_ASSERT(mid[i],>, i==0 ? min : top[i-1]);
		}
	}
	Axis(const double min, const double max, const unsigned nbins){
		this->min = min;
		top.resize(nbins);
		mid.resize(nbins);
		if(nbins==1) top[0] = max;
		else{
			double del = (max-min) / (double)nbins;
			for(int i=0; i<nbins; i++){
				top[i] = min + (i+1)*del;
				mid[i] = min + ((double)i + 0.5)*del;
			}
		}
	}

	Axis() {
		min = NaN;
	}

	unsigned size() const {return top.size();}

	int bin(const double x) const{
		if(x<min) return -1;
		else{
			// upper_bound returns first element greater than xval
			// values mark bin tops, so this is what we want
			int ind = upper_bound(top.begin(), top.end(), x) - top.begin();
			PRINT_ASSERT(ind,>=,0);
			PRINT_ASSERT(ind,<=,(int)top.size());
			return ind;
		}
	}

	double bottom(const int i) const{
		PRINT_ASSERT(i,<,size());
		return i==0 ? min : top[i-1];
	}

	double delta(const int i) const{
		PRINT_ASSERT(i,<,size());
		return top[i] - bottom(i);
	}
	double delta3(const int i) const{
		PRINT_ASSERT(i,<,size());
		return top[i]*top[i]*top[i] - bottom(i)*bottom(i)*bottom(i);
	}
	double max() const{
		return top[size()-1];
	}

	void write_HDF5(const string& name, H5::H5File file) const{
		hsize_t dims[1];
		H5::DataSpace dataspace;

		// write top
		dims[0] = size()+1;
		dataspace = H5::DataSpace(1,dims);
		H5::DataSet dataset = file.createDataSet(name+"[edge]",H5::PredType::IEEE_F64LE,dataspace);
		vector<double> tmp(size()+1);
		tmp[0] = min;
		for(unsigned i=1; i<size()+1; i++) tmp[i] = top[i-1];
		dataset.write(&tmp[0],H5::PredType::IEEE_F64LE);
		dataset.close();

		// write mid
		dims[0] = size();
		dataspace = H5::DataSpace(1,dims);
		dataset = file.createDataSet(name+"[mid]",H5::PredType::IEEE_F64LE,dataspace);
		dataset.write(&mid[0],H5::PredType::IEEE_F64LE);
		dataset.close();
	}
};

#endif
