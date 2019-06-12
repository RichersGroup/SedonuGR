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
#include "MomentSpectrumArray.h"

using namespace std;
namespace pc = physical_constants;


class TrajectoryData{
public:
	vector<double> ct, Ecom_Elab, Elab_Elab0, TMeV, Ye, rho;
	vector< vector< vector<double> > > Ndens, Fdens, Pdens;
	double nulab0;

	TrajectoryData(int NS, int NE){
		nulab0 = -1.e99;
		ct.resize(0);
		Ecom_Elab.resize(0);
		Elab_Elab0.resize(0);
		TMeV.resize(0);
		Ye.resize(0);
		rho.resize(0);

		Ndens.resize(NS);
		Fdens.resize(NS);
		Pdens.resize(NS);

		for(int s=0; s<NS; s++){
			Ndens[s].resize(NE);
			Fdens[s].resize(NE);
			Pdens[s].resize(NE);
			for(int g=0; g<NE; g++){
				Ndens[s][g].resize(0);
				Fdens[s][g].resize(0);
				Pdens[s][g].resize(0);
			}
		}
	}
};

void append_data(const Transport* sim, const EinsteinHelper* eh, double ct, TrajectoryData* td){
	double nulab = -eh->g.ndot(eh->kup);
	if(td->ct.size()==0) td->nulab0 = nulab;

	td->ct.push_back(ct);
	td->rho.push_back(sim->grid->rho.interpolate(eh->icube_vol));
	td->TMeV.push_back(sim->grid->T.interpolate(eh->icube_vol) * pc::k_MeV);
	td->Ye.push_back(sim->grid->Ye.interpolate(eh->icube_vol));
	td->Ecom_Elab.push_back(-eh->kup_tet[3]/eh->g.ndot(eh->kup));
	td->Elab_Elab0.push_back(nulab/td->nulab0);

	//MomentSpectrumArray<3>* dist = sim->grid->distribution[0];
	cout << "n=" << td->ct.size() << endl;
}

hid_t create_file(string filename, const TrajectoryData& td){
	hid_t file = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	hid_t file_space, dset, space1d;
	hsize_t ndims;

	// ct
	ndims = 1;
	hsize_t dims1[1] = {td.ct.size()};
	space1d = H5Screate_simple(ndims, dims1, dims1);

	dset = H5Dcreate(file, "ct(cm)", H5T_NATIVE_DOUBLE, space1d, H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
	H5Dwrite(dset, H5T_NATIVE_DOUBLE, space1d, space1d, H5P_DEFAULT, &td.ct[0]);

	dset = H5Dcreate(file, "rho(g|ccm)", H5T_NATIVE_DOUBLE, space1d, H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
	H5Dwrite(dset, H5T_NATIVE_DOUBLE, space1d, space1d, H5P_DEFAULT, &td.rho[0]);

	dset = H5Dcreate(file, "T(MeV)", H5T_NATIVE_DOUBLE, space1d, H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
	H5Dwrite(dset, H5T_NATIVE_DOUBLE, space1d, space1d, H5P_DEFAULT, &td.TMeV[0]);

	dset = H5Dcreate(file, "Ye", H5T_NATIVE_DOUBLE, space1d, H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
	H5Dwrite(dset, H5T_NATIVE_DOUBLE, space1d, space1d, H5P_DEFAULT, &td.Ye[0]);

	dset = H5Dcreate(file, "Ecom_Elab", H5T_NATIVE_DOUBLE, space1d, H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
	H5Dwrite(dset, H5T_NATIVE_DOUBLE, space1d, space1d, H5P_DEFAULT, &td.Ecom_Elab[0]);

	dset = H5Dcreate(file, "Elab_Elab0", H5T_NATIVE_DOUBLE, space1d, H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
	H5Dwrite(dset, H5T_NATIVE_DOUBLE, space1d, space1d, H5P_DEFAULT, &td.Elab_Elab0[0]);

	// clear resources
	H5Dclose(dset);
	H5Sclose(space1d);
	H5Fclose(file);
	return file;
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

	// open up the lua parameter file
	Lua lua;
	string script_file = string(argv[1]);
	lua.init( script_file );

	// set up the transport module (includes the grid)
	class testTransport : public Transport{
	public:
		void move(EinsteinHelper *eh, double* ct){
			PRINT_ASSERT(eh->ds_com,>=,0);
			PRINT_ASSERT(eh->N,>,0);
			PRINT_ASSERT(abs(eh->g.dot<4>(eh->kup,eh->kup)) / (eh->kup[3]*eh->kup[3]), <=, TINY);

			// save old values
			Tuple<double,4> old_kup = eh->kup;

			// convert ds_com into dlambda
			double dlambda = eh->ds_com / eh->kup_tet[3];
			PRINT_ASSERT(dlambda,>=,0);

			// get 2nd order x, 1st order estimate for k
			Tuple<double,4> order1 = old_kup * dlambda;
			for(size_t i=0; i<4; i++)
				eh->xup[i] += order1[i];
			if(DO_GR){
				Tuple<double,4> dk_dlambda = grid->dk_dlambda(*eh);
				Tuple<double,4> order2 = dk_dlambda * dlambda*dlambda * 0.5;
				eh->kup = old_kup + dk_dlambda * dlambda;
				for(size_t i=0; i<4; i++)
					eh->xup[i] += (abs(order2[i]/order1[i])<1. ? order2[i] : 0);
			}

			// get new background data
			update_eh_background(eh);
			if(eh->fate==moving)
				update_eh_k_opac(eh);

			double ds_com_new = dlambda*eh->kup_tet[3];
			*ct += (ds_com_new + eh->ds_com) / 2.;
		}
	};
	testTransport sim;
	sim.init(&lua);
	sim.reset_radiation();

	// read in the output data
	string recover_filename = lua.scalar<string>("recover_file");
	sim.grid->read_zones(recover_filename);
	const int NS = sim.grid->distribution.size();
	const int NE = sim.grid->nu_grid_axis.size();

	// read in starting points
	vector<vector<double> > xup(3), kup(3);
	xup[0] = lua.vector<double>("initial_x0");
	xup[1] = lua.vector<double>("initial_x1");
	xup[2] = lua.vector<double>("initial_x2");
	kup[0] = lua.vector<double>("initial_k0");
	kup[1] = lua.vector<double>("initial_k1");
	kup[2] = lua.vector<double>("initial_k2");
	int ntrajectories = xup[0].size();
	lua.close();


	for(int itraj=0; itraj<ntrajectories; itraj++){
		// initialize the EinsteinHelper
		EinsteinHelper eh;
		TrajectoryData td(NS,NE);
		for(int i=0; i<3; i++){
			eh.xup[i] = xup[i][itraj];
			eh.kup[i] = kup[i][itraj];
		}
		eh.xup[3] = 0;
		eh.s = 0;
		eh.N = 1;
		eh.N0 = eh.N;
		eh.fate = moving;

		sim.update_eh_background(&eh);
		eh.g.normalize_null_changeupt(eh.kup);
		sim.update_eh_k_opac(&eh);
		double ct = 0;
		append_data(&sim, &eh, ct, &td);
		while(eh.fate==moving){
			double d_zone = sim.grid->zone_min_length(eh.z_ind) / sqrt(Metric::dot_Minkowski<3>(eh.kup,eh.kup)) * eh.kup_tet[3];
			eh.ds_com = d_zone * sim.max_step_size;
			ct += eh.ds_com;
			sim.move(&eh, &ct);
			append_data(&sim, &eh, ct, &td);
		}
		create_file("trajectory"+to_string(itraj)+".h5", td);
	}

	// exit the program
	MPI_Finalize();
	return 0;
}
