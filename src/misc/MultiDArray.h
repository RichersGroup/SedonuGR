#ifndef _MULTIDARRAY_H
#define _MULTIDARRAY_H 1

#include <vector>
#include "global_options.h"
#include "Axis.h"
#include "mpi.h"

using namespace std;

//===================//
// InterpolationCube //
//===================//
template<unsigned ndims>
class InterpolationCube{
public:

	const static unsigned ncorners = (1<<ndims);
	constexpr static unsigned index(const unsigned LR[ndims]){
		unsigned result = 0;
		for(unsigned d=0; d<ndims; d++){
			PRINT_ASSERT(abs(LR[d]),<=,1);
			result = result << 1; // shift bit left
			result += LR[d]; // rightmost bit 1 if right value
		}
		return result;
	}
	constexpr static bool isRightIndex(unsigned i, unsigned d){ // 1 if right, 0 if left
		return bool( (i & ( 1 << d )) >> d );
	}
	double xLR[ndims][2];
	unsigned indices[ncorners];
	double weights[ncorners];
	double slope_weights[ndims][ncorners];

	InterpolationCube(){
		for(unsigned d=0; d<ndims; d++)	xLR[d][0] = xLR[d][1] = NaN;
		for(unsigned i=0; i<ncorners; i++){
			indices[i] = -1;
			weights[i] = NaN;
			for(unsigned d=0; d<ndims; d++)
				slope_weights[d][i] = NaN;
		}
	}

	bool inside_box(const double x[ndims]) const{
		bool result = true;
		for(unsigned d=0; d<ndims; d++){
			result = result && (x[d]>=xLR[d][0]);
			result = result && (x[d]<=xLR[d][1]);
		}
		return result;
	}

	void set_weights(const double x[ndims]){

		// total volume
		double V=1;
		for(unsigned d=0; d<ndims; d++){
			double dx = xLR[d][1] - xLR[d][0];
			V *= dx;
		}

		for(unsigned i=0; i<ncorners; i++){

			// volume/area associated with each point/line
			double dVol=1, dA[ndims];
			for(unsigned d=0; d<ndims; d++) dA[d] = 1.;
			for(unsigned d=0; d<ndims; d++){
				unsigned LR = not isRightIndex(i,d); // 1 if left, 0 if right
				double dx = x[d] - xLR[d][LR];
				dVol *= dx;
				for(unsigned d_deriv=0; d_deriv<ndims; d_deriv++)
					if(d_deriv != d) dA[d_deriv] *= dx;
			}

			// weights
			weights[i] = abs(dVol/V); // avoids if statement in above loop for sign
			for(unsigned d=0; d<ndims; d++){
				unsigned LR = isRightIndex(i,d);
				slope_weights[d][i] = (LR==1 ? 1.0 : -1.0) * abs(dA[d]/V);
			}
		}
	}
};


//=============//
// MultiDArray //
//=============//
template<typename T, unsigned nelements, unsigned int ndims>
class MultiDArray{
public:

	vector< Tuple<T,nelements> > y0;
	vector<Axis> axes;
	Tuple<unsigned int,ndims> stride;

	MultiDArray(){}

	void set_axes(const vector<Axis>& axes){
		this->axes = axes;
		PRINT_ASSERT(axes.size(),==,ndims);
		int size = 1;
		for(int i=ndims-1; i>=0; i--){
			stride[i] = size;
			size *= axes[i].size();
		}
		y0.resize(size);
		if(ndims==0) y0.resize(1);
	}

	MultiDArray<T,nelements,ndims> operator =(const MultiDArray<T,nelements,ndims>& input){
		PRINT_ASSERT(input.axes.size(),==,ndims);
		this->axes = input.axes;
		this->stride = input.stride;
		this->y0 = input.y0;
		return *this;
	}


	unsigned direct_index(const unsigned ind[ndims]) const{
		unsigned result = 0;
		for(unsigned i=0; i<ndims; i++){
			PRINT_ASSERT(ind[i],<,axes[i].size());
			result += ind[i]*stride[i];
		}
		PRINT_ASSERT(result,<,y0.size());
		return result;
	}
	void indices(const int z_ind, unsigned ind[ndims]) const{
		unsigned leftover=z_ind;
		PRINT_ASSERT(leftover,<,y0.size());
		for(unsigned i=0; i<ndims; i++){
			ind[i] = leftover / stride[i];
			leftover -= ind[i]*stride[i];
			PRINT_ASSERT(ind[i],<,axes[i].size());
		}
	}

	// get center value based on grid index
	const Tuple<double,nelements>& operator[](const unsigned i) const {
		return y0[i];
	}
	Tuple<T,nelements>& operator[](const unsigned i){
		return y0[i];
	}

	// dummy template allows it to compile with any value of NDIMS
	template<unsigned dummy>
	Tuple<T,nelements> interpolate(const InterpolationCube<dummy>& icube) const{
		PRINT_ASSERT(icube.ncorners,==,(1<<ndims));

		Tuple<T,nelements> result;
		result = y0[icube.indices[0]] * icube.weights[0];
		for(unsigned i=1; i<icube.ncorners; i++){
			PRINT_ASSERT(icube.indices[i],>=,0);
			Tuple<T,nelements> tmp = y0[icube.indices[i]] * icube.weights[i];
			result += tmp;
		}
		return result;
	}

	// dummy template allows it to compile with any value of NDIMS
	template<unsigned dummy>
	Tuple<Tuple<T,nelements>,ndims> interpolate_slopes(const InterpolationCube<dummy>& icube) const{
		PRINT_ASSERT(icube.ncorners,==,(1<<ndims));

		Tuple<Tuple<T,nelements>,ndims> result;
		for(unsigned d=0; d<ndims; d++){
			result[d] = y0[icube.indices[0]] * icube.slope_weights[d][0];;
			for(unsigned i=1; i<icube.ncorners; i++){
				PRINT_ASSERT(icube.indices[i],>=,0);
				result[d] += y0[icube.indices[i]] * icube.slope_weights[d][i];
			}
		}
		return result;
	}

	Tuple<T,nelements> slope(const unsigned z_ind, const unsigned direction) const{
		Tuple<T,nelements> result, yL, yR, y=y0[z_ind];
		yL = yR = NaN;
		double dxL=NaN, dxR=NaN;
		unsigned dir_ind[ndims];
		indices(z_ind,dir_ind);
		unsigned dir_indR[ndims], dir_indL[ndims];
		for(unsigned i=0; i<ndims; i++)
			dir_indR[i] = dir_indL[i] = dir_ind[i];
		dir_indR[direction] = dir_ind[direction]+1;
		dir_indL[direction] = dir_ind[direction]-1;

		if(dir_ind[direction] <= axes[direction].size()-2){
			yR = y0[direct_index(dir_indR)];                 // value at z_ind+1
			dxR = axes[direction].delta(dir_ind[direction]); // (top-bottom) of z_ind
		}
		if(dir_ind[direction] >= 1){
			yL = y0[direct_index(dir_indL)];                  // value at z_ind-1
			dxL = axes[direction].delta(dir_indL[direction]); // (top-bottom) of z_ind-1
		}

		PRINT_ASSERT(dxL>0,or,dxR>0);
		if(dxL>0 and dxR>0) result = ((y-yL)*dxR/dxL + (yR-y)*dxL/dxR)/(dxL+dxR);
		else if(dxL>0) result = (y-yL)/dxL;
		else result = (yR-y)/dxR;

		return result;
	}

	// dummy template allows it to compile with any value of NDIMS
	template<unsigned dummy>
	void set_InterpolationCube(InterpolationCube<dummy>* icube, const double x[ndims], const unsigned dir_ind_center[ndims]) const{
		if(not icube->inside_box(x) || ndims==0){

			// set boundary coordinates
			int dir_ind_left[ndims], dir_ind_right[ndims];
			for(unsigned d=0; d<ndims; d++){
				dir_ind_left[d] = (x[d]>axes[d].mid[dir_ind_center[d]] ? dir_ind_center[d] : dir_ind_center[d]-1);
				dir_ind_right[d] = dir_ind_left[d]+1;
				if(dir_ind_left[d] < 0){
					dir_ind_left[d] = 0;
					PRINT_ASSERT(dir_ind_right[d],==,dir_ind_left[d]);
					icube->xLR[d][0] = axes[d].min;
					icube->xLR[d][1] = axes[d].mid[dir_ind_right[d]];
				}
				else if(dir_ind_right[d] >= (int)axes[d].size()){
					dir_ind_right[d] = axes[d].size()-1;
					PRINT_ASSERT(dir_ind_left[d],==,dir_ind_right[d]);
					icube->xLR[d][0] = axes[d].mid[dir_ind_left[d]];
					icube->xLR[d][1] = axes[d].max();
				}
				else{
					icube->xLR[d][0] = axes[d].mid[dir_ind_left[d] ];
					icube->xLR[d][1] = axes[d].mid[dir_ind_right[d]];
				}

				// sanity checks
				PRINT_ASSERT(icube->xLR[d][1],>,icube->xLR[d][0]);
				PRINT_ASSERT(x[d],<=,icube->xLR[d][1]);
				PRINT_ASSERT(x[d],>=,icube->xLR[d][0]);
			}

			// set the global index of the values on the corners
			unsigned dir_ind[ndims];
			for(unsigned i=0; i<icube->ncorners; i++){
				for(unsigned d=0; d<ndims; d++)
					dir_ind[d] = (InterpolationCube<ndims>::isRightIndex(i,d) ? dir_ind_right[d] : dir_ind_left[d]);
				icube->indices[i] = direct_index(dir_ind);
			}
		}

		// calculate the weights associated with each corner
		icube->set_weights(x);
	}

	void wipe(){
		#pragma omp parallel for
		for(unsigned z=0; z<y0.size(); z++)
			y0[z] = 0;
	}

	unsigned size() const{
		return y0.size();
	}

	unsigned Ndims() const{
		return ndims;
	}

	void add(const unsigned ind[ndims], const Tuple<T,nelements> to_add){
		unsigned lin_ind = direct_index(ind);
		direct_add(lin_ind, to_add);
	}
	void direct_add(const unsigned lin_ind, const Tuple<T,nelements> to_add){
		for(unsigned i=0; i<nelements; i++){
			#pragma omp atomic
			y0[lin_ind][i] += to_add[i];
		}
	}

	void mpi_sum_scatter(vector<unsigned>& stop_list){
		PRINT_ASSERT(stop_list[stop_list.size()-1],==,y0.size());
		int MPI_nprocs, MPI_myID;
		MPI_Comm_size(MPI_COMM_WORLD, &MPI_nprocs);
		MPI_Comm_rank(MPI_COMM_WORLD, &MPI_myID);
		for(int p=0; p<MPI_nprocs; p++){
			const unsigned istart = (p==0 ? 0 : stop_list[p-1]);
			const unsigned ndoubles = (stop_list[p]-istart) * nelements;
			if(MPI_myID==p)
				MPI_Reduce(MPI_IN_PLACE, &y0[istart], ndoubles, MPI_DOUBLE, MPI_SUM, p, MPI_COMM_WORLD);
			else
				MPI_Reduce(&y0[istart],         NULL, ndoubles, MPI_DOUBLE, MPI_SUM, p, MPI_COMM_WORLD);
		}

		// and bring everything to proc0 as well
		mpi_gather(stop_list);
	}

	void mpi_sum(){
		int MPI_myID;
		MPI_Comm_rank(MPI_COMM_WORLD, &MPI_myID);
		if(MPI_myID==0)
			MPI_Reduce(MPI_IN_PLACE, &y0.front(), y0.size()*nelements, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		else
			MPI_Reduce(&y0.front(),         NULL, y0.size()*nelements, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	}

	void mpi_gather(vector<unsigned>& stop_list){
		int MPI_nprocs, MPI_myID;
		MPI_Comm_size(MPI_COMM_WORLD, &MPI_nprocs);
		MPI_Comm_rank(MPI_COMM_WORLD, &MPI_myID);
		PRINT_ASSERT((int)stop_list.size(),==,MPI_nprocs);
		PRINT_ASSERT(stop_list[stop_list.size()-1],==,y0.size());

		vector<int> sendcounts(MPI_nprocs), displs(MPI_nprocs);
		displs[0] = 0;
		sendcounts[0] = stop_list[0];
		for(int i=1; i<MPI_nprocs; i++){
			displs[i] = stop_list[i-1];
			sendcounts[i] = stop_list[i] - stop_list[i-1];
		}
		if(MPI_myID==0)
			MPI_Gatherv(MPI_IN_PLACE, -1, MPI_DOUBLE, &y0[0],&sendcounts.front(),&displs.front(), MPI_DOUBLE,0,MPI_COMM_WORLD);
		else
			MPI_Gatherv(&y0[displs[MPI_myID]], sendcounts[MPI_myID], MPI_DOUBLE, NULL, NULL, NULL, MPI_DOUBLE,0,MPI_COMM_WORLD);
	}

	void write_HDF5(H5::H5File file, const string name) {
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


//===================//
// ScalarMultiDArray //
//===================//
template<typename T, unsigned int ndims>
class ScalarMultiDArray : public MultiDArray<T,1,ndims>{
public:
	void add(const unsigned ind[ndims], const T to_add){
		unsigned lin_ind = this->direct_index(ind);
		direct_add(lin_ind, to_add);
	}
	void direct_add(const unsigned lin_ind, const T to_add){
		#pragma omp atomic
		this->y0[lin_ind][0] += to_add;
	}

	// get center value based on grid index
	const double& operator[](const unsigned i) const {
		return this->y0[i][0];
	}
	T& operator[](const unsigned i){
		return this->y0[i][0];
	}

	template<unsigned dummy>
	T interpolate(const InterpolationCube<dummy>& icube) const{
		return MultiDArray<T,1,ndims>::interpolate(icube)[0];
	}
	template<unsigned dummy>
	Tuple<double,ndims> interpolate_slopes(const InterpolationCube<dummy>& icube) const{
		Tuple<Tuple<T,1>,ndims> result;
		result = MultiDArray<T,1,ndims>::interpolate_slopes(icube);

		Tuple<T,ndims> return_value;
		for(unsigned i=0; i<ndims; i++) return_value[i] = result[i][0];
		return return_value;
	}
};

#endif
