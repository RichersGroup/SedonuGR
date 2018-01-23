#ifndef _MULTIDARRAY_H
#define _MULTIDARRAY_H 1

#include <vector>
#include "global_options.h"
#include "Axis.h"
#include "mpi.h"

using namespace std;

//=====================//
// INTERPOLATION ARRAY //
//=====================//
template<unsigned nelements, unsigned int ndims>
class MultiDArray{
public:

	vector< Tuple<double,nelements> > y0;
	vector<Axis> axes;
	vector< Tuple< Tuple<double,nelements> ,ndims> > dydx;
	Tuple<unsigned int,ndims> stride;

	MultiDArray(){}

	void set_axes(const vector<Axis>& axes){
		this->axes = axes;
		assert(axes.size()==ndims);
		int size = 1;
		for(int i=ndims-1; i>=0; i--){
			stride[i] = size;
			size *= axes[i].size();
		}
		y0.resize(size);
		if(ndims==0) y0.resize(1);
	}

	MultiDArray<nelements,ndims> operator =(const MultiDArray<nelements,ndims>& input){
		PRINT_ASSERT(input.axes.size(),==,ndims);
		this->axes = input.axes;
		this->stride = input.stride;
		this->dydx = input.dydx;
		this->y0 = input.y0;
		return *this;
	}


	unsigned direct_index(const unsigned ind[ndims]) const{
		int result = 0;
		for(int i=0; i<ndims; i++){
			PRINT_ASSERT(ind[i],<,axes[i].size());
			result += ind[i]*stride[i];
		}
		PRINT_ASSERT(result,<,y0.size());
		return result;
	}
	void indices(const int z_ind, unsigned ind[ndims]) const{
		unsigned leftover=z_ind;
		PRINT_ASSERT(leftover,<,y0.size());
		for(int i=0; i<ndims; i++){
			ind[i] = leftover / stride[i];
			leftover -= ind[i]*stride[i];
			PRINT_ASSERT(ind[i],<,axes[i].size());
		}
	}

	// get center value based on grid index
	const Tuple<double,nelements>& operator[](const unsigned i) const {return y0[i];}
	Tuple<double,nelements>& operator[](const unsigned i){return y0[i];}

	// get interpolated value
	Tuple<double,nelements> interpolate(const double x[ndims], const unsigned int ind[ndims]) const{
		unsigned z_ind = direct_index(ind);
		Tuple<double,nelements> result = y0[z_ind];
		if(dydx.size()>0) for(int i=0; i<ndims; i++){
			result += dydx[z_ind][i] * (x[i] - axes[i].mid[ind[i]]);
		}
		for(unsigned e=0; e<nelements; e++) PRINT_ASSERT(abs(result[e]),<,INFINITY);
		return result;
	}

	// set the slopes
	void calculate_slopes(const double minval, const double maxval){
		PRINT_ASSERT(maxval,>=,minval);
		dydx.resize(y0.size());

		#pragma omp parallel for
		for(unsigned int z=0; z<dydx.size(); z++){
			unsigned int ind[ndims], indp[ndims], indm[ndims];
			unsigned int zp, zm;
			double x, xp, xm;
			Tuple<double,nelements> y, yp, ym;
			double dxL, dxR;
			Tuple<double,nelements> dyL, dyR;
			Tuple<double,nelements> sL, sR;
			Tuple<double,nelements> slope;

			y=y0[z];
			indices(z, ind);
			for(unsigned int i=0; i<ndims; i++){
				if(axes[i].size()==1){
					dydx[z][i] = 0;
					continue;
				}

				// get the index for the plus and minus values
				for(unsigned int j=0; j<ndims; j++){
					indp[j] = ind[j];
					indm[j] = ind[j];
				}

				// get plus and minus values
				x=axes[i].mid[ind[i]];
				if(ind[i] > 0){
					indm[i] = ind[i]-1;
					zm = direct_index(indm);
					xm = axes[i].mid[indm[i]];
					ym = y0[zm];
					dxL = x-xm;
					dyL = y-ym;
					sL = dyL/dxL;
				}
				if(ind[i] < axes[i].size()-1){
					indp[i] = ind[i]+1;
					if(indp[i]>=axes[i].size()) sR=0;
					else{
						zp = direct_index(indp);
						xp = axes[i].mid[indp[i]];
						yp = y0[zp];
						dxR = xp-x;
						dyR = yp-y;
						sR = dyR/dxR;
					}
				}


				// get the actual slope
				if(ind[i]==0) slope = sR;
				else if(ind[i]==axes[i].size()-1) slope = sL;
				else slope = (sL*dxR + sR*dxL) / (dxR+dxL);

				dydx[z][i] = slope;
			}


			// check min/max values
			for(unsigned e=0; e<nelements; e++){
				double ybig = y[e];
				double ysmall = y[e];
				PRINT_ASSERT(y[e],<=,maxval);
				PRINT_ASSERT(y[e],>=,minval);

				for(unsigned i=0; i<ndims; i++){
					if(dydx[z][i][e] >= 0){
						ybig   += dydx[z][i][e] * (axes[i].top[ind[i]]    - axes[i].mid[ind[i]]);
						ysmall += dydx[z][i][e] * (axes[i].bottom(ind[i]) - axes[i].mid[ind[i]]);
					}
					else{
						ybig   += dydx[z][i][e] * (axes[i].bottom(ind[i]) - axes[i].mid[ind[i]]);
						ysmall += dydx[z][i][e] * (axes[i].top[ind[i]]    - axes[i].mid[ind[i]]);
					}
				}
				PRINT_ASSERT(ybig,>=,y[e]);
				PRINT_ASSERT(ysmall,<=,y[e]);

				for(unsigned i=0; i<ndims; i++){
					const double oldslope= dydx[z][i][e];
					if(ybig > maxval){
						dydx[z][i][e] = oldslope * (maxval-y[e]) / (ybig - y[e]);
						PRINT_ASSERT(ysmall,>=,minval);
					}
					if(ysmall < minval){
						dydx[z][i][e] = oldslope * (y[e]-minval) / (y[e]-ysmall);
						PRINT_ASSERT(ybig,<=,maxval);
					}
				}
			}
		}
	}

	void wipe(){
          #pragma omp parallel
	  {
	    #pragma omp for
	    for(unsigned z=0; z<y0.size(); z++)
	      y0[z] = 0;
            #pragma omp for collapse(2)
	    for(unsigned z=0; z<dydx.size(); z++)
	      for(unsigned i=0; i<ndims; i++)
		dydx[z][i] = 0;
	  }
	}

	unsigned size() const{
		return y0.size();
	}

	unsigned Ndims() const{
		return ndims;
	}

	void add(const unsigned ind[ndims], const Tuple<double,nelements> to_add){
		unsigned lin_ind = direct_index(ind);
		direct_add(lin_ind, to_add);
	}
	void direct_add(const unsigned lin_ind, const Tuple<double,nelements> to_add){
		for(unsigned i=0; i<nelements; i++){
			#pragma omp atomic
			y0[lin_ind][i] += to_add[i];
		}
	}

	void MPI_combine()
	{
		int MPI_myID;
		MPI_Comm_rank( MPI_COMM_WORLD, &MPI_myID);
		MPI_Request request;
		const int tag = 0;

		// average the flux (receive goes out of scope after section)
		if(MPI_myID==0) MPI_Reduce(MPI_IN_PLACE, &y0.front(), size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		else            MPI_Reduce(&y0.front(),  &y0.front(), size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	}

	void MPI_AllCombine()
	{
		int MPI_myID;
		MPI_Comm_rank( MPI_COMM_WORLD, &MPI_myID);
		MPI_Request request;
		const int tag = 0;

		// average the flux (receive goes out of scope after section)
		MPI_Allreduce(MPI_IN_PLACE, &y0.front(), size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	}

	void write_HDF5(H5::H5File file, const string name) const {
		hsize_t dims[ndims+1];
		for(unsigned i=0; i<ndims; i++) dims[i] = axes[i].size(); // number of bins
		H5::DataSpace dataspace;
		if(nelements==1)
			dataspace = H5::DataSpace(ndims,dims);
		else{
			dims[ndims] = nelements;
			dataspace = H5::DataSpace(ndims+1,dims);
		}
		H5::DataSet dataset = file.createDataSet(name,H5::PredType::IEEE_F64LE,dataspace);

		// write the data (converting to single precision)
		// assumes phi increases fastest, then mu, then nu
		dataset.write(&y0.front(), H5::PredType::IEEE_F64LE);
		dataset.close();
	}
};

template<unsigned int ndims>
class ScalarMultiDArray : public MultiDArray<1,ndims>{
public:
	void add(const unsigned ind[ndims], const double to_add){
		unsigned lin_ind = this->direct_index(ind);
		direct_add(lin_ind, to_add);
	}
	void direct_add(const unsigned lin_ind, const double to_add){
		#pragma omp atomic
		this->y0[lin_ind][0] += to_add;
	}

	// get center value based on grid index
	const double& operator[](const unsigned i) const {return this->y0[i][0];}
	double& operator[](const unsigned i){return this->y0[i][0];}

	// get interpolated value
	double interpolate(const double x[ndims], const unsigned int ind[ndims]) const{
		Tuple<double,1> result = MultiDArray<1,ndims>::interpolate(x,ind);
		return result[0];
	}
};

#endif
