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
#include "Lua.h"
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
	const int rank0 = (MPI_myID == 0);

	// open up the lua parameter file
	Lua lua;
	string script_file = string(argv[1]);
	lua.init( script_file );
	double r_sch = lua.scalar<double>("Grid1DSchwarzschild_r_sch");

	class testTransport : public Transport{
	public:
		void set_particle(EinsteinHelper& eh){
			particles.resize(1);
			particles[0] = eh.p;
		}
		virtual void move(EinsteinHelper *eh){
			for(unsigned i=0; i<4; i++) cout << eh->p.xup[i] << "\t";
			for(unsigned i=0; i<4; i++) cout << eh->p.kup[i] << "\t";
			for(unsigned i=0; i<4; i++) cout << eh->kup_tet[i] << "\t";
			cout << eh->nu() << endl;

			Transport::move(eh);
			if(eh->p.xup[0]<0 and eh->p.xup[1]<0) eh->p.fate = absorbed;
		}
	};

	// set up the transport module (includes the grid)
	testTransport sim;
	sim.init(&lua);

	// start only one neutrino
	EinsteinHelper eh;

	vector<double> xup,kup;
	xup.resize(4);
	kup.resize(4);
	xup = lua.vector<double>("initial_xup");
	kup = lua.vector<double>("initial_kup");
	eh.p.xup[3] = 0;
	for(int i=0; i<4; i++){
		eh.p.xup[i] = xup[i];
		eh.p.kup[i] = kup[i];
	}
	eh.p.s = 0;
	eh.p.N = 1;
	eh.p.fate = moving;

	sim.update_eh_background(&eh);
	eh.g.normalize_null_preserveupt(eh.p.kup);
	sim.update_eh_k_opac(&eh);

	sim.grid->interpolate_opacity(&eh);

	sim.set_particle(eh);
	sim.step();
	sim.write(0);

	// read in time stepping parameters
	lua.close();

	// exit the program
	MPI_Finalize();
	return 0;
}
