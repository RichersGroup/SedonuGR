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
#include "global_options.h"
#include "physical_constants.h"
#include "Species.h"
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
}

void Species::init(Lua* lua, Transport* simulation)
{
	// initialize MPI parallelism
	int my_rank;
	MPI_Comm_rank( MPI_COMM_WORLD, &my_rank );
	const int rank0 = (my_rank == 0);

	// set the pointer to see the simulation info
	sim = simulation;
	PRINT_ASSERT(sim->grid->z.size(),>,0);

	// call child's init function
	myInit(lua);

	// allocate space for the grid eas spectrum containers
	if(rank0) cout << "#   Setting up the eas arrays...";
	abs_opac.resize(sim->grid->z.size());
	scat_opac.resize(sim->grid->z.size());
	emis.resize(sim->grid->z.size());
	biased_emis.resize(sim->grid->z.size());

	// allocate space for each eas spectrum
	PRINT_ASSERT(nu_grid.size(),>,0);
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

    // set up the spectrum in each zone
	if(rank0) cout << "#   Setting up the distribution function...";
	int n_mu = lua->scalar<int>("distribution_nmu");
	int n_phi = lua->scalar<int>("distribution_nphi");

	// temporary spectrum to be used for distribution function initialization
	LocateArray tmp_mugrid, tmp_phigrid;
	if(n_mu>=0)	tmp_mugrid.init(-1,1,n_mu);
	else{
		vector<double> bin_tops = lua->vector<double>("distribution_mugrid_interfaces");
		bin_tops.push_back(1);
		tmp_mugrid.init(bin_tops,-1);
	}
	if(n_mu>=0)	tmp_phigrid.init(-pc::pi, pc::pi, n_phi);
	else{
		vector<double> bin_tops = lua->vector<double>("distribution_phigrid_interfaces");
		bin_tops.push_back(pc::pi);
		tmp_phigrid.init(bin_tops,-pc::pi);
	}

	// use the above to initialize the zone's distribution function
	SpectrumArray tmp_spectrum;
	tmp_spectrum.init(nu_grid, tmp_mugrid, tmp_phigrid);

	//#pragma omp parallel for
	for(unsigned z_ind=0; z_ind<sim->grid->z.size(); z_ind++){
		sim->grid->z[z_ind].distribution.push_back(tmp_spectrum);
	}
	if(rank0) cout << "finished." << endl;

	// set nugrid interpolation method
	int imeth  = lua->scalar<int>("opac_interp_method");
	if     (imeth == 0) nu_grid.interpolation_method = constant;
	else if(imeth == 1) nu_grid.interpolation_method = linear;
	else if(imeth == 2) nu_grid.interpolation_method = logarithmic;
	else				nu_grid.interpolation_method = power;
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
	double imp = interpolate_importance(nu,zone_index);
	PRINT_ASSERT(imp,>,0);
	*Ep *= 1.0/imp;
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

double Species::bin_emis(const int z_ind, const int g) const{
	return emis[z_ind].get_value(g) * emis[z_ind].N;
}

unsigned Species::number_of_bins(){
	return nu_grid.size();
}

double Species::sum_opacity(const int z_ind, const int group) const{
	return abs_opac[z_ind][group] + scat_opac[z_ind][group];
}
