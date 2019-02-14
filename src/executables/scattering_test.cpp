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
		void testgrid(EinsteinHelper *eh){
			for(size_t igin=0; igin<grid->nu_grid_axis.size(); igin++){
				for(size_t igout=0; igout<grid->nu_grid_axis.size(); igout++){
					grid->scat_opac[0][igin]=1.0;
					grid->scat_opac[1][igin]=1.0;
					grid->scat_opac[2][igin]=1.0;
					grid->scattering_delta[0][igin]=0.0;
					grid->scattering_delta[1][igin]=0.0;
					grid->scattering_delta[2][igin]=0.0;
					if(igin==igout){
					grid->partial_scat_opac[0][igin][igin]=1.0;
					grid->partial_scat_opac[1][igin][igin]=1.0;
					grid->partial_scat_opac[2][igin][igin]=1.0;
					} else {
					grid->partial_scat_opac[0][igin][igout]=0.0;
					grid->partial_scat_opac[1][igin][igout]=0.0;
					grid->partial_scat_opac[2][igin][igout]=0.0;
					}
				}
			}
		}
	};

	// set up the transport module (includes the grid)
	testScatter sim;
	sim.init(&lua);

	// start 1 neutrino
	EinsteinHelper eh;
	ofstream myfile;
	eh.kup[0]=1*pc::h;
	eh.kup[1]=0;
	eh.kup[2]=0;
	eh.kup[3]=sim.grid->nu_grid_axis.mid[0]*pc::h;
	eh.xup[0]=0;
	eh.xup[1]=0;
	eh.xup[2]=0;
	eh.xup[3]=0;
	//eh.dir_ind[0]=1;
	eh.s = 0;
	eh.N = 1;
	eh.N0 = eh.N;
	//cout<<eh.kup_tet<<endl;
	sim.testgrid(&eh);
	cout<<"set testgrid"<<endl;
	sim.update_eh_background(&eh);
	cout<<"calling update ehkopac"<<endl;
        sim.update_eh_k_opac(&eh); 
	cout<<"updated ehkopac"<<endl;
	sim.scatter(&eh);
	cout<<"back from scatter"<<endl;
  	myfile.open ("elastic_isotropic_kernel.dat");
	myfile<<eh.kup[0]<<" "<<eh.kup[1]<<" "<<eh.kup[2]<<" "<<eh.kup[3]<<"\n";		
	myfile.close();
	
	// read in time stepping parameters
	lua.close();

	// exit the program
	MPI_Finalize();
	return 0;
}
