#ifndef _SPECTRUM_ARRAY_H
#define _SPECTRUM_ARRAY_H 1

#include <string>
#include "global_options.h"
#include "particle.h"
#include "locate_array.h"
#include "transport.h"

using namespace std;

// default values
#define DEFAULT_NAME "spectrum_array"

class spectrum_array {

private:

	// bin arrays
	// values represent bin upper walls (the single locate_array.min value is the leftmost wall)
	// underflow is combined into leftmost bin (right of the locate_array.min)
	// overflow is combined into the rightmost bin (left of locate_array[size-1])
	locate_array time_grid;
	locate_array wave_grid;
	locate_array mu_grid;
	locate_array phi_grid;

	// counting arrays
	vector<double> flux;
	vector<int>    click;

	// Indexing
	int index(const int t,const int l,const int m,const int p) const;

public:

	// Initialize
	void init(const locate_array tg, const locate_array wg, const locate_array mg, const locate_array pg);
	void init(const std::vector<double>,const std::vector<double>,const int,const int);

	// MPI functions
	void MPI_average();

	// Count a packets
	void count(const double t, const double w, const double E, const vector<double> D);

	//  void normalize();
	void rescale(const double);
	void wipe();

	// Print out
	void print(const int iw, const int species) const;
};

#endif
