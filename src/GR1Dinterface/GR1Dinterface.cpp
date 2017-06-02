#include "Transport.h"
#include "GridGR1D.h"
#include "Neutrino_GR1D.h"
#include "mpi.h"
#include <iostream>

using namespace std;

extern "C"
void initialize_gr1d_sedonu_(const double *x1i, const int* n_GR1D_zones, const int* M1_imaxradii, const int* ghosts1, Transport** sim){

	// initialize MPI parallelism
	int my_rank,n_procs;
	MPI_Init(NULL,NULL);
	MPI_Comm_rank( MPI_COMM_WORLD, &my_rank );
	MPI_Comm_size( MPI_COMM_WORLD, &n_procs);
	const int rank0 = (my_rank == 0);


	int c_M1_imaxradii = *M1_imaxradii-1;
	cout << "#   M1_imaxradii = " << c_M1_imaxradii << endl;
	cout << "#   ghosts1 = " << *ghosts1 << endl;
	cout << "#   n_GR1D_zones = " << *n_GR1D_zones << endl;
	cout << "#   radius range: " << x1i[*ghosts1] << " - " << x1i[c_M1_imaxradii+1] << endl;
	cout << "#   next lowest radius: " << x1i[c_M1_imaxradii] << endl;

	// declare the transport class and save it to the fortran module
	*sim = new Transport;

	// open up the lua parameter file
	Lua lua("param.lua");

	// set up the grid
	const int nzones = *M1_imaxradii - *ghosts1;
	GridGR1D* grid = new GridGR1D;
	grid->initialize_grid(x1i,nzones,*ghosts1);
	(*sim)->grid = grid;

	// set up the transport module (includes the grid)
	(*sim)->init(&lua);
	for(int s=0; s<(*sim)->species_list.size(); s++){
		Neutrino_GR1D* tmpSpecies =static_cast<Neutrino_GR1D*>((*sim)->species_list[s]);
		tmpSpecies->ghosts1 = *ghosts1;
		tmpSpecies->n_GR1D_zones = *n_GR1D_zones;
	}
}

extern "C"
void calculate_mc_closure_(double* q_M1, double* q_M1_extra, const double* eas,
		const double* rho, const double* T, const double* Ye, const double* v1,
		const int* iter, Transport** sim){

	// set velocities and opacities
	static_cast<GridGR1D*>((*sim)->grid)->set_fluid(rho, T, Ye, v1);

	for(int s=0; s<(*sim)->species_list.size(); s++)
		static_cast<Neutrino_GR1D*>((*sim)->species_list[s])->set_eas_external(eas);

	// do MC calculation
	(*sim)->step();
}
