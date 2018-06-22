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
#include "LuaRead.h"
#include "Transport.h"

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
	int MPI_myID;
	MPI_Init( &argc, &argv );
	MPI_Comm_rank( MPI_COMM_WORLD, &MPI_myID );
	const int rank0 = (MPI_myID == 0);

	// start timer
	double proc_time_start = MPI_Wtime();

	// open up the lua parameter file
	Lua lua;
	string script_file = ( argc>1 ? string(argv[1]) : "param.lua");
	lua.init( script_file );

	// set up the transport module (includes the grid)
	Transport sim;
	sim.init(&lua);

	// read in time stepping parameters
	int max_n_iter  = lua.scalar<int>("max_n_iter");\
	if(max_n_iter < 0){
		max_n_iter = MAXLIM;
	}
	double max_time_hours = lua.scalar<double>("max_time_hours");
	double max_time_seconds = max_time_hours * 3600.0;
	if(max_time_seconds < 0){
		max_time_seconds = INFINITY;
	}
	lua.close();

	// initial output
	sim.write(0);

	//===========//
	// TIME LOOP //
	//===========//
	for(int it=1; it<=max_n_iter; it++)
	{
		double time_now = MPI_Wtime();
		if(rank0) cout << "# Elapsed Time: " << (time_now - proc_time_start) / 60. << " minutes" << endl;
		if(time_now - proc_time_start <= max_time_seconds){
			// do transport step
			sim.step();

			// printout time step
			sim.write(it);
		}
		else break;
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
