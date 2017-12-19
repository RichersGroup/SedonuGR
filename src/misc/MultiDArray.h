#include <vector>
#include "global_options.h"
#include "Axis.h"

using namespace std;

class MultiDInterface{
public:
	virtual ~MultiDInterface() {}
	virtual double get(const unsigned int ind[]) const=0;
	virtual void set(const unsigned int ind[], const double y)=0;
	virtual double interpolate(const double x[], const unsigned int ind[]) const=0;
	virtual void calculate_slopes()=0;
};


//=====================//
// INTERPOLATION ARRAY //
//=====================//
template<unsigned int ndims>
class MultiDArray : public MultiDInterface{
private:

	vector<double> y0;
	vector<Tuple<double,ndims> > dydx;
	Tuple<unsigned int,ndims> stride;
	Tuple<Axis*,ndims> axes;

	unsigned int index(const unsigned int ind[ndims]) const{
		int result = 0;
		for(int i=0; i<ndims; i++) result += ind[i]*stride[i];
		return result;
	}
	void indices(const int z_ind, unsigned int ind[ndims]) const{
		for(int i=0; i<ndims; i++) ind[i] = z_ind % stride[i];
	}

public:

	MultiDArray(const vector<Axis*>& axes){
		assert(axes.size()==ndims);
		int size = 1;
		for(int i=0; i<ndims; i++){
			this->axes[i] = axes[i];
			stride[ndims-1-i] = size;
			size *= axes[i]->size();
		}
		y0.resize(size);
	}


	// get center value based on grid index
	double get(const unsigned int ind[ndims]) const{
		return y0[index(ind)];
	}
	void set(const unsigned int ind[ndims], const double y){
		y0[index(ind)] = y;
	}

	// get interpolated value
	double interpolate(const double x[ndims], const unsigned int ind[ndims]) const{
		unsigned int z_ind = index(ind);
		double result = y0[z_ind];
		if(dydx.size()>0) for(int i=0; i<ndims; i++)
			result += dydx[z_ind][i] * (x[i] - axes[i]->mid[ind[i]]);
		return result;
	}

	// set the slopes
	void calculate_slopes(){
		unsigned int ind[ndims], indp[ndims], indm[ndims];
		unsigned int zp, zm;
		double x, xp, xm;
		double y, yp, ym;
		double dxL, dxR;
		double dyL, dyR;
		double sL, sR;
		double slope;

		dydx.resize(y0.size());
		for(unsigned int z=0; z<dydx.size(); z++){
			indices(z, ind);
			for(unsigned int i=0; i<ndims; i++){

				// get the index for the plus and minus values
				for(unsigned int j=0; j<ndims; j++){
					indp[j] = ind[j];
					indm[j] = ind[j];
				}
				indp[i]++;
				indm[i]--;
				zp = index(indp);
				zm = index(indm);

				// get plus and minus values
				y=y0[z];
				x=axes[i]->mid[ind[i]];
				if(ind[i] > 0){
					xm = axes[i]->mid[indm[i]];
					ym = y0[zm];
					dxL = x-xm;
					dyL = y-ym;
					sL = dyL/dxL;
				}


				// get the actual slope
				if(ind[i]==0) slope = sR;
				else if(ind[i]==axes[i]->size()-1) slope = sL;
				else slope = (dxR*sL + dxL*sR) / (dxR+dxL);
				dydx[z][i] = slope;
			}
		}
	}
};
