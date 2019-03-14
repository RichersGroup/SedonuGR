#include "MultiDArray.h"
#include <iostream>

using namespace std;

bool print_test(const double result, const double expected){
	cout << result << " ";
	bool pass = fabs(result-expected) / (fabs(expected)>0 ? fabs(expected) : 1.0) < TINY;
	if(pass) cout << endl;
	else cout << "\tFAIL: expected " << expected << endl;
	return pass;
}

int main(){
	bool pass = true;

	const int nbins=4, ndims=3;
	cout << "Axis i bottom mid top" << endl;
	Axis xAxis(-4,4,nbins);
	for(int i=0; i<nbins; i++){
		cout << i << " " << xAxis.bottom(i) << " " << xAxis.mid[i] << " " << xAxis.top[i] << endl;
	}
	cout << "mda values: " << endl;
	vector<Axis> axes(ndims,xAxis);
	MultiDArray<double,ndims,ndims> mda; // 3 dims and 3 components on each dim
	mda.set_axes(axes);
	for(int i=0; i<mda.size(); i++){
		cout << i << ": {";
		for(int j=0; j<ndims; j++){
			mda[i][j] = i * (j+1);
			cout << mda[i][j] << " ";
		}
		cout << "}" << endl;
	}

	cout << "|=============|" << endl;
	cout << "| INDEX TESTS |" << endl;
	cout << "|=============|" << endl;
	size_t ind[ndims] = {1,2,3};
	cout << "input: { ";
	for(int i=0; i<ndims; i++) cout << ind[i] << " ";
	cout << "}" << endl;
	cout << "direct index: ";
	size_t direct_index = mda.direct_index(ind);
	pass = pass and print_test(direct_index, 27);
	size_t out_ind[ndims];
	mda.indices(direct_index, out_ind);
	for(int i=0; i<ndims; i++){
		cout << "index" << i<<": ";
		pass = pass and print_test(out_ind[i], ind[i]);
	}

	cout << "|====================|" << endl;
	cout << "| INTERPOLATION_CUBE |" << endl;
	cout << "|====================|" << endl;
	size_t dir_ind_center[ndims] = {3,2,2}; //{1,1,1}
	double x[ndims] = {3.5,1+TINY,1.5}; //{0,0,0}
	cout << "input: {";
	for(int i=0; i<ndims; i++) cout << x[i] << " ";
	cout << "}" << endl;
	size_t expected_indices[8] = {58,58, 62,62, 59,59, 63,63};//{21,37, 25,41, 22,38, 26,42};
	vector<double> expected_weights{3./8,3./8, 0,0, 1./8,1./8, 0,0};
	vector<vector<double> > expected_slope_weights{
		{-3./4 , 3./4 ,      0,     0, -1./4 , 1./4 ,      0,     0},
		{-3./16,-3./16,  3./16, 3./16, -1./16,-1./16,  1./16, 1./16},
		{-1./4 ,-1./4 ,      0,     0,  1./4 , 1./4 ,      0,     0},
	};

	InterpolationCube<ndims> icube;
	mda.set_InterpolationCube(&icube, x, dir_ind_center);
	icube.set_slope_weights(x);
	cout << "Inside cube: ";
	pass = pass and print_test(icube.inside_box(x), true);
	cout << "icube indices:" << endl;
	for(int i=0; i<icube.ncorners; i++){
		pass = pass and print_test(icube.indices[i], expected_indices[i]);
	}
	cout << "interpolate: "<<endl;
	Tuple<double,ndims> interpval = mda.interpolate(icube);
	pass = pass and print_test(interpval[0], 58.25);
	pass = pass and print_test(interpval[1], 116.5);
	pass = pass and print_test(interpval[2], 174.75);
	cout << "icube weights:" << endl;
	for(int i=0; i<icube.ncorners; i++){
		pass = print_test(icube.weights[i], expected_weights[i]) and pass;
	}
	Tuple<Tuple<double,ndims>,ndims> interpslope = mda.interpolate_slopes(icube);
	cout << "interpolate slope:"<<endl;
	pass = pass and print_test(interpslope[0][0], 0);
	pass = pass and print_test(interpslope[0][1], 0);
	pass = pass and print_test(interpslope[0][2], 0);
	pass = pass and print_test(interpslope[1][0], 2);
	pass = pass and print_test(interpslope[1][1], 4);
	pass = pass and print_test(interpslope[1][2], 6);
	pass = pass and print_test(interpslope[2][0], .5);
	pass = pass and print_test(interpslope[2][1], 1);
	pass = pass and print_test(interpslope[2][2], 1.5);
	for(int d=0; d<ndims; d++){
		cout << "icube slope weights " << d <<":"<<endl;
		for(int i=0; i<icube.ncorners; i++){
			pass = print_test(icube.slope_weights[d][i], expected_slope_weights[d][i]) and pass;
		}
	}

	assert(pass);
	return 0;
}
