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
#include "nulib_interface.h"
#include "Species.h"
#include "Grid.h"
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
	class testScatter : public Transport{
	public:
		void testgrid(){
		  for(size_t s=0; s<species_list.size(); s++){
		  			for(size_t igin=0; igin<10; igin++){
		  			grid->scat_opac[s][igin] = 0;
					grid->abs_opac[s][igin] = 0;
		  			for(size_t igout=0; igout<10; igout++){
		  				grid->scattering_delta[s][igin]=0.0;
		  				grid->partial_scat_opac[s][igin][igout] = (igin==igout ? 1 : 0);
		  				grid->scat_opac[s][igin] += grid->partial_scat_opac[s][igin][igout];
		  			}
		  		}
		  	}
		}
	};

	// set up the transport module (includes the grid)
	testScatter sim;
	sim.init(&lua);
	lua.close();

	// start 1 neutrino
	EinsteinHelper eh;
	eh.kup[0]=0.33*sim.grid->nu_grid_axis.mid[0]*pc::h;
	eh.kup[1]=0;
	eh.kup[2]=0;
	eh.kup[3]=0.33*sim.grid->nu_grid_axis.mid[0]*pc::h;
	eh.xup[0]=0;
	eh.xup[1]=0;
	eh.xup[2]=0;
	eh.xup[3]=0;
	eh.s = 0;
	eh.N = 1;
	eh.N0 = eh.N;
	sim.testgrid();
	sim.update_eh_background(&eh);
        sim.update_eh_k_opac(&eh); 

	// run many copies of the same neutrino and output the results.
	ofstream myfile;
  	myfile.open ("elastic_isotropic_kernel.dat");
	for(int i=0; i<1000; i++){
	  EinsteinHelper tmp = eh;
	  sim.scatter(&tmp);
	  myfile<<tmp.kup[0]<<" "<<tmp.kup[1]<<" "<<tmp.kup[2]<<" "<<tmp.kup[3]<<"\n";
	}
	myfile.close();
	

	// exit the program
	MPI_Finalize();
	return 0;
}
