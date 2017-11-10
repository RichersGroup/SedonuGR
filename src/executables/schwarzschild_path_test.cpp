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

double run_test(const int nsteps, const bool rank0, const double rho, const double T_MeV, const double target_ye, Transport& sim, ofstream& outf){
	if(rank0) cout << endl << "Currently running: rho=" << rho << "g/ccm T_core=" << T_MeV << "MeV Ye=" << target_ye << endl;

	// set the fluid properties
	sim.grid->z[0].rho = rho;
	double error = 1.2;
	if(sim.equilibrium_T>0 ) sim.grid->z[0].T  = min(T_MeV/pc::k_MeV*error,sim.T_max);
	else               sim.grid->z[0].T  = T_MeV/pc::k_MeV;
	if(sim.equilibrium_Ye>0) sim.grid->z[0].Ye = min(target_ye*error,sim.Ye_max);
	else               sim.grid->z[0].Ye = target_ye;

	sim.grid->z[0].T = T_MeV/pc::k_MeV; //min(T_MeV/pc::k_MeV*1.1,100/pc::k_MeV);
	sim.grid->z[0].Ye = target_ye; //min(target_ye*1.1,0.55);
	double T_core = T_MeV/pc::k_MeV;

	// reconfigure the core
	double munue = nulib_eos_munue(rho,T_core,target_ye);
	sim.init_core(sim.r_core,T_core,munue);
	PRINT_ASSERT(sim.core_species_luminosity.N,>,0);

	// check max optical depth
	double max_opac = 0;
	for(unsigned z_ind=0; z_ind<sim.grid->z.size(); z_ind++){
		for(unsigned s=0; s<sim.species_list.size(); s++){
			for(unsigned g=0; g<sim.species_list[s]->number_of_bins(); g++){
				double opac = sim.species_list[s]->sum_opacity(z_ind,g);
				if(opac>max_opac) max_opac = opac;
			}
		}
	}
	double optical_depth = max_opac * sim.grid->zone_min_length(0);
	if(rank0) cout << " Optical Depth: " << optical_depth << endl;

	// do the transport step
	for(int i=0; i<nsteps; i++) sim.step();

	// write the data out to file
	outf << rho << "\t" << T_MeV << "\t" << target_ye << "\t" << munue*pc::ergs_to_MeV << "\t";
	if(rank0) sim.grid->write_line(outf,0);
	return optical_depth;
}

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
	string script_file = string("param.lua");
	lua.init( script_file );
	double r_sch = lua.scalar<double>("Grid1DSchwarzschild_r_sch");

	class testTransport : public Transport{
	public:
		void set_particle(Particle p){
			particles.resize(1);
			particles[0] = p;
			particles[0].print();
		}
	};

	// set up the transport module (includes the grid)
	testTransport sim;
	sim.init(&lua);

	// start only one neutrino
	Particle p;
	double r0 = 1.5*r_sch;
	p.s = 0;
	p.e = 1;
	p.xup[0] = r0;
	p.xup[1] = M_PI/2.0;
	p.xup[2] = 0;
	p.xup[3] = 0;
	p.kup[0] = (1.0 - r_sch/r0);
	p.kup[1] = 0;
	p.kup[2] = 0;
	p.kup[3] = 1;
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
