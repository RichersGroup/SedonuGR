/*
//  Copyright (c) 2015, California Institute of Technology and the Regents
//  of the University of California, based on research sponsored by the
//  United States Department of Energy. All rights reserved.
//
//  This file is part of Sedonu.
//
//  Sedonu is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  Neither the name of the California Institute of Technology (Caltech)
//  nor the University of California nor the names of its contributors 
//  may be used to endorse or promote products derived from this software
//  without specific prior written permission.
//
//  Sedonu is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with Sedonu.  If not, see <http://www.gnu.org/licenses/>.
//
*/

#include <mpi.h>
#include <string>
#include "Lua.h"
#include "transport.h"

using namespace std;

//--------------------------------------------------------
// The main code
// The user writes this for their own needs
//--------------------------------------------------------
int main(int argc, char **argv)
{
	//============//
	// INITIALIZE //
	//============//
	// initialize MPI parallelism
	int my_rank,n_procs;
	MPI_Init( &argc, &argv );
	MPI_Comm_rank( MPI_COMM_WORLD, &my_rank );
	MPI_Comm_size( MPI_COMM_WORLD, &n_procs);
	const int rank0 = (my_rank == 0);

	// start timer
	double proc_time_start = MPI_Wtime();

	// open up the lua parameter file
	Lua lua;
	string script_file = ( argc>1 ? string(argv[1]) : "param.lua");
	lua.init( script_file );

	// set up the transport module (includes the grid)
	transport sim;
	sim.init(&lua);

	// read in time stepping parameters
	int max_n_iter  = lua.scalar<int>("max_n_iter");
	lua.close();

	// initial output
	sim.write(0);

	//===========//
	// TIME LOOP //
	//===========//
	if (rank0) printf("%12s %12s %12s %12s\n","iteration","t","dt","n_particles");
	for(int it=1; it<=max_n_iter; it++)
	{
		// do transport step
		sim.step();

		// printout time step
		sim.write(it);
		if(rank0) printf("%12d %12d\n",it,sim.total_particles());
	}

	//===================//
	// FINALIZE AND EXIT //
	//===================//
	// calculate the elapsed time
	double proc_time_end = MPI_Wtime();
	double time_wasted = proc_time_end - proc_time_start;
	if (rank0) printf("#\n# CALCULATION took %.3e seconds or %.3f mins or %.3f hours\n",
			time_wasted,time_wasted/60.0,time_wasted/60.0/60.0);

	// exit the program
	MPI_Finalize();
	return 0;
}
