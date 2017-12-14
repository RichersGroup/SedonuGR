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
	T_min = NaN;
	T_max = NaN;
	Ye_min = NaN;
	Ye_max = NaN;
	rho_min = NaN;
	rho_max = NaN;
	sim = NULL;
	ID = MAXLIM;
	cutoff=0;
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

	// basic parameters
	cutoff = lua->scalar<double>("cdf_cutoff");

	//===================================//
    // set neutrino's min and max values //
	//===================================//
    T_min   =  0.0;
    T_max   =  numeric_limits<double>::infinity();
    Ye_min  = 0;
    Ye_max  = 1;
    rho_min = -numeric_limits<double>::infinity();
    rho_max =  numeric_limits<double>::infinity();

	//============================//
	// CALL CHILD'S INIT FUNCTION //
	//============================//
	myInit(lua);

    //==========================//
    // setup for frequency grid // only if child didn't
    //==========================//
    double minval = 0;
    double trash, tmp=0;
    vector<double> bintops = vector<double>(0);

    // get the frequency grid
    if(nu_grid.size()==0){
    	int n_nu   = lua->scalar<int>("nugrid_n");
    	if(n_nu>0){
    		double nu_start = lua->scalar<double>("nugrid_start");
    		double nu_stop  = lua->scalar<double>("nugrid_stop");
    		PRINT_ASSERT(nu_stop,>,nu_start);
    		nu_grid.init(nu_start/pc::h_MeV,nu_stop/pc::h_MeV,n_nu);
    	}
    	else{
    		string nugrid_filename = lua->scalar<string>("nugrid_filename");
    		ifstream nugrid_file;
    		nugrid_file.open(nugrid_filename.c_str());
    		nugrid_file >> trash >> minval;
    		minval /= pc::h_MeV;
    		while(nugrid_file >> trash >> tmp){
    			tmp /= pc::h_MeV;
    			if(bintops.size()>0) PRINT_ASSERT(tmp,>,bintops[bintops.size()-1]);
    			else PRINT_ASSERT(tmp,>,minval);
    			bintops.push_back(tmp);
    		}
    		nugrid_file.close();
            nu_grid.init(bintops,minval,nu_grid.interpolation_method);
    	}
	}
	PRINT_ASSERT(nu_grid.size(),>,0);

    // set the interpolation method
    int imeth  = lua->scalar<int>("opac_interp_method");
	if     (imeth == 0) nu_grid.interpolation_method = constant;
	else if(imeth == 1) nu_grid.interpolation_method = linear;
	else if(imeth == 2) nu_grid.interpolation_method = logarithmic;
	else				nu_grid.interpolation_method = power;

    //===========================//
	// intialize output spectrum // only if child didn't
	//===========================//
	LocateArray tmp_mugrid, tmp_phigrid;
	if(spectrum.size()==0){
		int nmu  = lua->scalar<int>("spec_n_mu");
		int nphi = lua->scalar<int>("spec_n_phi");
		tmp_mugrid.init( -1     , 1     , nmu );
		tmp_phigrid.init(-pc::pi, pc::pi, nphi);
		spectrum.init(nu_grid, tmp_mugrid, tmp_phigrid);
	}
	PRINT_ASSERT(spectrum.size(),>,0);

	//============================================//
	// allocate space for the grid eas containers //
	//============================================//
	if(rank0) cout << "#   Setting up the eas arrays...";
	abs_opac.resize(sim->grid->z.size());
	scat_opac.resize(sim->grid->z.size());
	emis.resize(sim->grid->z.size());
	biased_emis.resize(sim->grid->z.size());
	if(sim->use_scattering_kernels==1){
		normalized_phi0.resize(sim->grid->z.size());
		scattering_delta.resize(sim->grid->z.size());
		for(int z_ind=0; z_ind<sim->grid->z.size(); z_ind++){
			normalized_phi0[z_ind].resize(nu_grid.size());
			scattering_delta[z_ind].resize(nu_grid.size());
			for(int igin=0; igin<nu_grid.size(); igin++){
				normalized_phi0[z_ind][igin].resize(nu_grid.size());
				scattering_delta[z_ind][igin].resize(nu_grid.size(),0);
			}
		}
	}

	// allocate space for each eas spectrum
	if(sim->n_emit_core>0 || sim->n_emit_core_per_bin>0) core_emis.resize(nu_grid.size());
	int iorder = lua->scalar<int>("cdf_interpolation_order");
    //#pragma omp parallel for
	for(unsigned i=0; i<abs_opac.size();  i++){
		abs_opac[i].resize(nu_grid.size());
		scat_opac[i].resize(nu_grid.size());
		emis[i].resize(nu_grid.size());
		emis[i].interpolation_order = iorder;
		biased_emis[i].resize(nu_grid.size());
		biased_emis[i].interpolation_order = iorder;
	}
	if(rank0) cout << "finished." << endl;

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
	      ((PolarSpectrumArray*)tmp_spectrum)->init(nu_grid, tmp_mugrid, tmp_phigrid);
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
	    ((PolarSpectrumArray*)tmp_spectrum)->init(nu_grid, tmp_mugrid, tmp_phigrid);
	    PRINT_ASSERT(((PolarSpectrumArray*)tmp_spectrum)->size(),>,0);
	    
	  }
	  
	  //-- MOMENT SPECTRUM --------------------
	  else if(distribution_type == "Moments"){
	    int order = lua->scalar<int>("distribution_moment_order");
	    tmp_spectrum = new MomentSpectrumArray;
	    ((MomentSpectrumArray*)tmp_spectrum)->init(nu_grid, order);
	  }

	  //-- RADIAL MOMENT SPECTRUM --------------------
	  else if(distribution_type == "RadialMoments"){
	    int order = lua->scalar<int>("distribution_moment_order");
	    tmp_spectrum = new RadialMomentSpectrumArray;
	    ((RadialMomentSpectrumArray*)tmp_spectrum)->init(nu_grid, order);
	  }

	  //-- RADIAL MOMENT SPECTRUM --------------------
	  else if(distribution_type == "GR1D"){
	    tmp_spectrum = new GR1DSpectrumArray;
	    ((GR1DSpectrumArray*)tmp_spectrum)->init(nu_grid);
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


//-----------------------------------------------------
// set cdf to blackbody distribution
// units of emis.N: erg/s/cm^2/ster
//-----------------------------------------------------
void Species::set_cdf_to_BB(const double T, const double chempot, CDFArray& emis){
    #pragma omp parallel for ordered
	for(unsigned j=0;j<nu_grid.size();j++)
	{
		double nu  = nu_grid.center(j);
		double dnu = nu_grid.delta(j);
        #pragma omp ordered
		emis.set_value(j, blackbody(T,chempot,nu)*dnu);
	}
	emis.normalize();
}


//----------------------------------------------------------------
// return a randomly sampled frequency
// for a particle emitted from the core or zone
//----------------------------------------------------------------
double Species::interpolate_importance(double nu, const int z_ind) const{
	PRINT_ASSERT(z_ind,>=,0);
	PRINT_ASSERT(z_ind,<,(int)emis.size());
	PRINT_ASSERT(nu,>=,nu_grid.min);
	PRINT_ASSERT(nu,<=,nu_grid[nu_grid.size()-1]);

	// frequency-integrated biasing already taken care of in emit_zones by changing average particle
	// energy depending on the luminosity of the zone and the number of particles emitted from it.
	double result = biased_emis[z_ind].interpolate_pdf(nu,&nu_grid)/**biased_emis[z_ind].N*/ / (emis[z_ind].interpolate_pdf(nu,&nu_grid)/**emis[z_ind].N*/);
	return result;
}
double Species::sample_core_nu(const int g) const
{
	PRINT_ASSERT(nu_grid.min,>=,0);
	return sample_nu(core_emis);
}
double Species::sample_zone_nu(const int zone_index, double *Ep, const int g) const
{
	PRINT_ASSERT(nu_grid.min,>=,0);
	double nu = sample_nu(biased_emis[zone_index],g);
	if(sim->importance_bias>0){
		double imp = interpolate_importance(nu,zone_index);
		PRINT_ASSERT(imp,>,0);
		*Ep *= 1.0/imp;
	}
	PRINT_ASSERT(nu,>,0);
	return nu;
}
double Species::sample_nu(const CDFArray& input_emis, const int g) const{
	// randomly pick a frequency
	double rand = sim->rangen.uniform();

	// make sure we don't get a zero frequency
	if(rand==0 && nu_grid.min==0) while(rand==0) rand = sim->rangen.uniform();

	return input_emis.invert(rand,&nu_grid,g);
}

//----------------------------------------------------------------
// return the emissivity integrated over nu for the core (erg/s)
//----------------------------------------------------------------
double Species::integrate_core_emis() const
{
	return core_emis.N;
}

//----------------------------------------------------------------
// return the emissivity integrated over nu for a zone (erg/s/ster/cm^3)
//----------------------------------------------------------------
double Species::integrate_zone_emis(const int zone_index) const
{
	return emis[zone_index].N;
}
double Species::integrate_zone_biased_emis(const int zone_index) const{
	return biased_emis[zone_index].N;
}


//----------------------------------------------------------------
// return the lepton emissivity integrated over nu for a zone (#/s/ster/cm^3)
// ASSUMES linear cdf sampling
//----------------------------------------------------------------
double Species::integrate_zone_lepton_emis(const int zone_index) const
{
	double l_emis = 0;
	for(unsigned i=0; i<emis[zone_index].size(); i++)
	{
		l_emis += lepton_number * emis[zone_index].get_value(i) / (pc::h*nu_grid.x[i]);
	}
	return l_emis * emis[zone_index].N;
}

//-----------------------------------------------------------------
// get opacity at the frequency
//-----------------------------------------------------------------
void Species::get_opacity(const double com_nu, const int z_ind, double* a, double* s) const
{
	PRINT_ASSERT(z_ind,>=,-1);
	PRINT_ASSERT(com_nu,>,0);

	// absorption and scattering opacities
	*a = max(nu_grid.value_at(com_nu, abs_opac[z_ind]),0.0);
	*s = max(nu_grid.value_at(com_nu,scat_opac[z_ind]),0.0);

	PRINT_ASSERT(*a,>=,0);
	PRINT_ASSERT(*s,>=,0);
	PRINT_ASSERT(*a,<,INFINITY);
	PRINT_ASSERT(*s,<,INFINITY);
}


double Species::bin_emis(const int z_ind, const int g) const{
	return emis[z_ind].get_value(g) * emis[z_ind].N;
}

unsigned Species::number_of_bins(){
	return nu_grid.size();
}

double Species::sum_opacity(const int z_ind, const int group) const{
	return abs_opac[z_ind][group] + scat_opac[z_ind][group];
}

//-------------------------------------------------------------
// Sample outgoing neutrino direction and energy
//-------------------------------------------------------------
void Species::sample_scattering_final_state(const int z_ind, LorentzHelper &lh, const double cosTheta) const{
	PRINT_ASSERT(sim->use_scattering_kernels,>,0);
	PRINT_ASSERT(scattering_delta.size(),>,0);
	PRINT_ASSERT(normalized_phi0.size(),>,0);

	// get ingoing energy index
	int igin = nu_grid.locate(lh.p_nu());
	if(igin == nu_grid.size()) igin--;

	// get outgoing energy
	double out_nu = normalized_phi0[z_ind][igin].invert(sim->rangen.uniform(),&nu_grid);
	lh.scale_p_energy(out_nu/lh.p_nu());

	// get outgoing energy index
	int igout = nu_grid.locate(out_nu);

	// bias outgoing direction to be isotropic. Very inefficient for large values of delta.
	double delta = scattering_delta[z_ind][igin][igout];
	PRINT_ASSERT(fabs(delta),<,3.0);
	if(fabs(delta)<=1.0) lh.scale_p_number(1.0 + delta*cosTheta);
	else{
		double b = 2.*fabs(delta) / (3.-fabs(delta));
		if(delta>1.0) lh.scale_p_number( pow(1.+cosTheta, b) );
		else          lh.scale_p_number( pow(1.-cosTheta, b) );
	}


	// sample outgoing direction
	//double U = sim->rangen.uniform();
	//double tiny = 1.e-4;
	//if(abs(delta)<1.e-4) // isotropic
	//	*outCosTheta = 2.*U - 1.;
	//else{
	//	double b2m4ac = sqrt(1. - delta*(2.-delta-4.*U));
	//	double sgn = (delta<0 ? -1.0 : 1.0);
	//	*outCosTheta = 1./delta * (-1. + sgn*b2m4ac);
	//}
}
