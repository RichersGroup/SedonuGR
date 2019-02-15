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
#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include "LuaRead.h"
#include "Transport.h"
#include "Species.h"
#include "Grid.h"
#include "nulib_interface.h"
#include "global_options.h"

using namespace std;
namespace pc = physical_constants;

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

	// open up the lua parameter file
	Lua lua;
	string script_file = string(argv[1]);
	lua.init( script_file );

	class testTransport : public Transport{
	public:
		void move(EinsteinHelper *eh){
			Transport::move(eh);
			if(eh->xup[0]<0 or eh->xup[1]<0) eh->fate = absorbed;

			for(size_t i=0; i<4; i++) cout << eh->xup[i] << "\t";
			for(size_t i=0; i<4; i++) cout << eh->kup[i] << "\t";
			for(size_t i=0; i<4; i++) cout << eh->kup_tet[i] << "\t";
			cout << eh->g.gtt << "\t";
			cout << eh->nu() << endl;
		}
	};

	// set up the transport module (includes the grid)
	testTransport sim;
	sim.init(&lua);
	sim.reset_radiation();

	// start only one neutrino
	EinsteinHelper eh;

	vector<double> xup,kup;
	xup.resize(4);
	kup.resize(4);
	xup = lua.vector<double>("initial_xup");
	kup = lua.vector<double>("initial_kup");
	eh.xup[3] = 0;
	for(int i=0; i<4; i++){
		eh.xup[i] = xup[i];
		eh.kup[i] = kup[i];
	}
	eh.s = 0;
	eh.N = 1;
	eh.N0 = eh.N;
	eh.fate = moving;

	sim.update_eh_background(&eh);
	sim.update_eh_k_opac(&eh);
	ParticleEvent event;
	while(eh.fate==moving){
		sim.which_event(&eh,&event);
		sim.move(&eh);
	}

	// read in time stepping parameters
	lua.close();

	// exit the program
	MPI_Finalize();
	return 0;
}
