#include "Transport.h"
#include "GridGR1D.h"
#include "mpi.h"
#include <iostream>

using namespace std;

extern "C"
void initialize_gr1d_sedonu_(const double *x1i, const int* M1_imaxradii, const int* ghosts1,
		const double* nulibtable_ebottom, const double* nulibtable_etop, const int* nulibtable_number_groups,
		const int* nulibtable_number_species, Transport* sim){

	// initialize MPI parallelism
	int my_rank,n_procs;
	MPI_Init(NULL,NULL);
	MPI_Comm_rank( MPI_COMM_WORLD, &my_rank );
	MPI_Comm_size( MPI_COMM_WORLD, &n_procs);
	const int rank0 = (my_rank == 0);


	int c_M1_imaxradii = *M1_imaxradii-1;
	cout << "M1_imaxradii = " << c_M1_imaxradii << endl;
	cout << "ghosts1 = " << *ghosts1 << endl;
	cout << "nulibtable_number_groups = " << *nulibtable_number_groups << endl;
	cout << "number_species_to_evolve = " << *nulibtable_number_species << endl;
	cout << "radius range: " << x1i[*ghosts1] << " - " << x1i[c_M1_imaxradii+1] << endl;
	cout << "next lowest radius: " << x1i[c_M1_imaxradii] << endl;

	// declare the transport class and save it to the fortran module
	sim = new Transport;

	// open up the lua parameter file
	cout << "=== WORKING ON LUA ===" << endl;
	Lua lua("param.lua");

	// set up the grid
	cout << "=== WORKING ON GRID ===" << endl;
	const int nzones = *M1_imaxradii - *ghosts1;
	GridGR1D* grid = new GridGR1D;
	grid->initialize_grid(&(x1i[*ghosts1+1]),nzones);
	sim->grid = grid;

	// set up the transport module (includes the grid)
	cout << "=== WORKING ON TRANSPORT ===" << endl;
	sim->init(&lua);

	exit(0);
}

extern "C"
void calculate_mc_closure_(double* q_M1, double* q_M1_extra, const double* eas, const double* v1, const double* vphi1, const int* iter, Transport* sim){

	// set velocities and opacities

	// do MC calculation
	sim->step();

	std::cout << "ASDFASDFASDFASDFASDFASDF" << std::endl;
	exit(0);
}
