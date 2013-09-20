#ifndef _SPECTRUM_ARRAY_H
#define _SPECTRUM_ARRAY_H 1

#include <string>
#include "particle.h"
#include "locate_array.h"

using namespace std;

// default values
#define DEFAULT_NAME "spectrum_array"

class spectrum_array {
 
private:

  // spectrum_array name
  char name[1000];

  // bin arrays
  // values represent bin walls
  // underflow is combined into leftmost bin (right of the leftmost wall)
  // overflow is stored separately right of rightmost wall
  locate_array time_grid;
  locate_array wave_grid;
  locate_array mu_grid;
  locate_array phi_grid;

  // counting arrays
  std::vector<double> flux;
  std::vector<int>    click;

  // Indexing
  int n_elements;
  int a1, a2, a3;
  int index(int,int,int,int);
    
public:

  // constructors
  spectrum_array();
  
  // Initialize
  void init(std::vector<double>,std::vector<double>,std::vector<double>,std::vector<double>);
  void init(std::vector<double>,std::vector<double>,int,int);
  void log_init(std::vector<double>,std::vector<double>,int,int);
  void set_name(const char *n);

  // MPI functions
  void MPI_average();
  
  // Count a packets
  void count(double t, double w, double E, double *D);

  //  void normalize();
  void rescale(double);
  void wipe();

  // Print out
  void print();
};

#endif
