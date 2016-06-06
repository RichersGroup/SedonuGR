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
	int my_rank,n_procs;
	MPI_Init( &argc, &argv );
	MPI_Comm_rank( MPI_COMM_WORLD, &my_rank );
	MPI_Comm_size( MPI_COMM_WORLD, &n_procs);
	const int rank0 = (my_rank == 0);

	// read command line input
	double max_logrho, min_logrho, rho0;
	double max_logT  , min_logT  , T0;
	double max_ye    , min_ye    , ye0;
	int n_rho, n_T, n_ye;
	PRINT_ASSERT(argc,==,15);
	sscanf(argv[2 ], "%lf", &min_logrho);
	sscanf(argv[3 ], "%lf", &max_logrho);
	sscanf(argv[4 ], "%lf", &rho0);
	sscanf(argv[5 ], "%d" , &n_rho);
	sscanf(argv[6 ], "%lf", &min_logT);
	sscanf(argv[7 ], "%lf", &max_logT);
	sscanf(argv[8 ], "%lf", &T0);
	sscanf(argv[9 ], "%d" , &n_T);
	sscanf(argv[10], "%lf", &min_ye);
	sscanf(argv[11], "%lf", &max_ye);
	sscanf(argv[12], "%lf", &ye0);
	sscanf(argv[13], "%d" , &n_ye);
	double dlogT   = (max_logT   - min_logT  ) / ((double)n_T   - 1.0);
	double dlogrho = (max_logrho - min_logrho) / ((double)n_rho - 1.0);
	double dye     = (max_ye     - min_ye    ) / ((double)n_ye  - 1.0);

	// start timer
	double proc_time_start = MPI_Wtime();

	// read in eos table
	nulib_eos_read_table(argv[14]);

	// open up the lua parameter file
	Lua lua;
	string script_file = string(argv[1]);
	lua.init( script_file );

	// set up the transport module (includes the grid)
	Transport sim;
	sim.init(&lua);

	// read in time stepping parameters
	int max_n_iter  = lua.scalar<int>("max_n_iter");
	lua.close();

	// open the output file
	ofstream outf;
	if(rank0) outf.open("results.dat");

	//==================//
	// TEMPERATURE LOOP //
	//==================//
	double max_optical_depth = 0;
	for(int i=0; i<n_T; i++){
		double logT = min_logT + i*dlogT;
		double optical_depth = run_test(max_n_iter, rank0, rho0,pow(10,logT),ye0,sim,outf);
		max_optical_depth = (optical_depth>max_optical_depth ? optical_depth : max_optical_depth);
	}
	//=========//
	// YE LOOP //
	//=========//
	for(int i=0; i<n_ye; i++){
		double ye = min_ye + i*dye;
		double optical_depth = run_test(max_n_iter, rank0,rho0,T0,ye,sim,outf);
		max_optical_depth = (optical_depth>max_optical_depth ? optical_depth : max_optical_depth);
	}
	//==============//
	// DENSITY LOOP //
	//==============//
	for(int i=0; i<n_rho; i++){
		double logrho = min_logrho + i*dlogrho;
		double optical_depth = run_test(max_n_iter, rank0,pow(10,logrho),T0,ye0,sim,outf);
		max_optical_depth = (optical_depth>max_optical_depth ? optical_depth : max_optical_depth);
	}

	//===================//
	// FINALIZE AND EXIT //
	//===================//
	// calculate the elapsed time
	if(rank0) cout << "MAXIMUM OPTICAL DEPTH: " << max_optical_depth << endl;
	double proc_time_end = MPI_Wtime();
	double time_wasted = proc_time_end - proc_time_start;
	if (rank0) printf("#\n# CALCULATION took %.3e seconds or %.3f mins or %.3f hours\n",
			time_wasted,time_wasted/60.0,time_wasted/60.0/60.0);

	// exit the program
	if(rank0) outf.close();
	MPI_Finalize();
	return 0;
}
