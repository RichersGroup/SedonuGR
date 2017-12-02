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
		void set_particle(Particle p){
			particles.resize(1);
			particles[0] = p;
		}
		void move(LorentzHelper *lh, int *z_ind){
			lh->particle_readonly(lab)->print();
			Transport::move(lh,z_ind);
			if(lh->p_xup()[0]<0 and lh->p_xup()[1]<0) lh->set_p_fate(absorbed);
		}
	};

	// set up the transport module (includes the grid)
	testTransport sim;
	sim.init(&lua);

	// start only one neutrino
	Particle p;

	vector<double> xup,kup;
	xup.resize(4);
	kup.resize(4);
	xup = lua.vector<double>("initial_xup");
	kup = lua.vector<double>("initial_kup");
	for(int i=0; i<4; i++){
		p.xup[i] = xup[i];
		p.kup[i] = kup[i];
	}

	p.s = 0;
	p.e = 1;
	sim.grid->normalize_null(p.kup,4,p.xup);
	p.tau = INFINITY;
	p.fate = moving;
	sim.set_particle(p);
	sim.step();
	sim.write(0);

	// read in time stepping parameters
	lua.close();

	// open the output file
	ofstream outf;
	if(rank0) outf.open("results.dat");

	// exit the program
	if(rank0) outf.close();
	MPI_Finalize();
	return 0;
}
