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
#include <fstream>
#include "global_options.h"
#include "physical_constants.h"
#include "Species.h"
#include "PolarSpectrumArray.h"
#include "MomentSpectrumArray.h"
#include "RadialMomentSpectrumArray.h"
#include "GR1DSpectrumArray.h"
#include "Transport.h"
#include "Grid.h"

using namespace std;
namespace pc = physical_constants;

Species::Species(){
	weight = NaN;
	lepton_number = MAXLIM;
	sim = NULL;
	ID = MAXLIM;
	nu_grid=NULL;
}

void Species::init(Lua* lua, Transport* simulation)
{
	// initialize MPI parallelism
	int MPI_myID;
	MPI_Comm_rank( MPI_COMM_WORLD, &MPI_myID );
	const int rank0 = (MPI_myID == 0);

	// set the pointer to see the simulation info
	sim = simulation;
	PRINT_ASSERT(sim->grid->z.size(),>,0);
	nu_grid = &(sim->grid->nu_grid);


	//============================//
	// CALL CHILD'S INIT FUNCTION //
	//============================//
	myInit(lua);

    // set the interpolation method
    int imeth  = lua->scalar<int>("opac_interp_method");
	if     (imeth == 0) nu_grid->interpolation_method = constant;
	else if(imeth == 1) nu_grid->interpolation_method = linear;
	else if(imeth == 2) nu_grid->interpolation_method = logarithmic;
	else				nu_grid->interpolation_method = power;

    //===========================//
	// intialize output spectrum // only if child didn't
	//===========================//
	LocateArray tmp_mugrid, tmp_phigrid;
	if(spectrum.size()==0){
		int nmu  = lua->scalar<int>("spec_n_mu");
		int nphi = lua->scalar<int>("spec_n_phi");
		tmp_mugrid.init( -1     , 1     , nmu );
		tmp_phigrid.init(-pc::pi, pc::pi, nphi);
		spectrum.init(*nu_grid, tmp_mugrid, tmp_phigrid);
	}
	PRINT_ASSERT(spectrum.size(),>,0);

	//============================================//
	// allocate space for the grid eas containers //
	//============================================//
	if(rank0) cout << "#   Setting up the eas arrays...";
	abs_opac.resize(sim->grid->z.size());
	scat_opac.resize(sim->grid->z.size());
	emis.resize(sim->grid->z.size());
	if(sim->use_scattering_kernels==1){
		normalized_phi0.resize(sim->grid->z.size());
		scattering_delta.resize(sim->grid->z.size());
		for(int z_ind=0; z_ind<sim->grid->z.size(); z_ind++){
			normalized_phi0[z_ind].resize(nu_grid->size());
			scattering_delta[z_ind].resize(nu_grid->size());
			for(int igin=0; igin<nu_grid->size(); igin++){
				normalized_phi0[z_ind][igin].resize(nu_grid->size());
				scattering_delta[z_ind][igin].resize(nu_grid->size(),0);
			}
		}
	}

	// allocate space for each eas spectrum
	int iorder = lua->scalar<int>("cdf_interpolation_order");
    //#pragma omp parallel for
	for(unsigned i=0; i<abs_opac.size();  i++){
		abs_opac[i].resize(nu_grid->size());
		scat_opac[i].resize(nu_grid->size());
		emis[i].resize(nu_grid->size());
		emis[i].interpolation_order = iorder;
	}
	if(rank0) cout << "finished." << endl;

    double minval = 0;
    double trash, tmp=0;
    vector<double> bintops = vector<double>(0);

    //==================================//
    // set up the spectrum in each zone //
    //==================================//
	if(rank0) cout << "#   Setting up the distribution function...";
	string distribution_type = lua->scalar<string>("distribution_type");
	for(unsigned z_ind=0; z_ind<sim->grid->z.size(); z_ind++){
	  SpectrumArray* tmp_spectrum;
	  
	  //-- POLAR SPECTRUM -----------------
	  if(distribution_type == "Polar"){
	    tmp_spectrum = new PolarSpectrumArray;

	    int n_mu = lua->scalar<int>("distribution_nmu");
	    int n_phi = lua->scalar<int>("distribution_nphi");
	    if(n_mu>0 && n_phi>0){
	      tmp_mugrid.init(-1,1,n_mu);
	      tmp_phigrid.init(-pc::pi, pc::pi, n_phi);
	      ((PolarSpectrumArray*)tmp_spectrum)->init(*nu_grid, tmp_mugrid, tmp_phigrid);
	    }
	    else{
	      // mu ---------------
	      string mugrid_filename = lua->scalar<string>("distribution_mugrid_filename");
	      ifstream mugrid_file;
	      mugrid_file.open(mugrid_filename.c_str());
	      mugrid_file >> trash >> minval;
	      bintops = vector<double>(0);
	      while(mugrid_file >> trash >> tmp){
		if(bintops.size()>0) PRINT_ASSERT(tmp,>,bintops[bintops.size()-1]);
		else PRINT_ASSERT(tmp,>,minval);
		bintops.push_back(tmp);
	      }
	      bintops[bintops.size()-1] = 1.0;
	      PRINT_ASSERT(minval,==,-1);
	      tmp_mugrid.init(bintops,minval);
	      
	      // phi --------------
	      string phigrid_filename = lua->scalar<string>("distribution_phigrid_filename");
	      ifstream phigrid_file;
	      phigrid_file.open(phigrid_filename.c_str());
	      phigrid_file >> trash >> minval;
	      minval -= pc::pi;
	      bintops = vector<double>(0);
	      while(phigrid_file >> trash >> tmp){
		tmp -= pc::pi;
		if(bintops.size()>0) PRINT_ASSERT(tmp,>,bintops[bintops.size()-1]);
		else PRINT_ASSERT(tmp,>,minval);
		bintops.push_back(tmp);
	      }
	      bintops[bintops.size()-1] = pc::pi;
	      PRINT_ASSERT(minval,==,-pc::pi);
	      tmp_phigrid.init(bintops,minval);
	    }
	    ((PolarSpectrumArray*)tmp_spectrum)->init(*nu_grid, tmp_mugrid, tmp_phigrid);
	    PRINT_ASSERT(((PolarSpectrumArray*)tmp_spectrum)->size(),>,0);
	    
	  }
	  
	  //-- MOMENT SPECTRUM --------------------
	  else if(distribution_type == "Moments"){
	    int order = lua->scalar<int>("distribution_moment_order");
	    tmp_spectrum = new MomentSpectrumArray;
	    ((MomentSpectrumArray*)tmp_spectrum)->init(*nu_grid, order);
	  }

	  //-- RADIAL MOMENT SPECTRUM --------------------
	  else if(distribution_type == "RadialMoments"){
	    int order = lua->scalar<int>("distribution_moment_order");
	    tmp_spectrum = new RadialMomentSpectrumArray;
	    ((RadialMomentSpectrumArray*)tmp_spectrum)->init(*nu_grid, order);
	  }

	  //-- RADIAL MOMENT SPECTRUM --------------------
	  else if(distribution_type == "GR1D"){
	    tmp_spectrum = new GR1DSpectrumArray;
	    ((GR1DSpectrumArray*)tmp_spectrum)->init(*nu_grid);
	  }

	  //-- CATCH ------------------
	  else{
	    cout << "Distribution type not found." << endl;
	    assert(0);
	  }
	  int distribution_polar_basis = lua->scalar<int>("distribution_polar_basis");
	  tmp_spectrum->rotated_basis = distribution_polar_basis;
	  
	  // push a distribution function to each zone
	  sim->grid->z[z_ind].distribution.push_back(tmp_spectrum);
	  sim->grid->z[z_ind].Edens_com.push_back(0);
	  sim->grid->z[z_ind].Ndens_com.push_back(0);
	  PRINT_ASSERT(sim->grid->z[z_ind].distribution.size(),==,ID+1);
	}
	
	if(rank0) cout << "finished." << endl;
}


//-----------------------------------------------------------------
// get opacity at the frequency
//-----------------------------------------------------------------
void Species::get_opacity(const double com_nu, const int z_ind, double* a, double* s) const
{
	PRINT_ASSERT(z_ind,>=,-1);
	PRINT_ASSERT(com_nu,>,0);

	// absorption and scattering opacities
	*a = max(nu_grid->value_at(com_nu, abs_opac[z_ind]),0.0);
	*s = max(nu_grid->value_at(com_nu,scat_opac[z_ind]),0.0);

	PRINT_ASSERT(*a,>=,0);
	PRINT_ASSERT(*s,>=,0);
	PRINT_ASSERT(*a,<,INFINITY);
	PRINT_ASSERT(*s,<,INFINITY);
}

