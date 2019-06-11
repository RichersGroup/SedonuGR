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
#include <iostream>
#include <iomanip>
#include <fstream>
#include "Grid.h"
#include "LuaRead.h"
#include "nulib_interface.h"
#include "global_options.h"
#include "Transport.h"
#include "H5Cpp.h"
#include "Neutrino_GR1D.h"
#include "MomentSpectrumArray.h"
#include "RadialMomentSpectrumArray.h"
#include "GR1DSpectrumArray.h"

using namespace std;
namespace pc = physical_constants;

Grid::Grid(){
	xAxes.resize(NDIMS);
	sim = NULL;
	do_annihilation=0;
	tetrad_rotation = cartesian;
}
//------------------------------------------------------------
// initialize the grid
//------------------------------------------------------------
void Grid::init(Lua* lua, Transport* insim)
{
	// MPI stuff
	int MPI_myID;
	MPI_Comm_rank( MPI_COMM_WORLD, &MPI_myID );
	bool rank0 = (MPI_myID==0);

	// set the transport pointer
	sim = insim;

	// read the model file or fill in custom model
	read_model_file(lua);
	for(size_t i=0; i<rho.size(); i++){
		// modify problematic fluid quantities
		if(rho[i] > sim->rho_max){
			if(sim->verbose) cout << "WARNING: resetting rho["<<i<<"] from "<<rho[i]<<" to "<<sim->rho_max<<endl;
			rho[i] = sim->rho_max;
		}
		// if(rho[i] < sim->rho_min){
		// 	if(sim->verbose) cout << "WARNING: resetting rho["<<i<<"] from "<<rho[i]<<" to "<<sim->rho_min<<endl;
		// 	rho[i] = sim->rho_min;
		// }
		if(T[i] > sim->T_max){
			if(sim->verbose) cout << "WARNING: resetting T["<<i<<"] from "<<T[i]<<" to "<<sim->T_max<<endl;
			T[i] = sim->T_max;
		}
		if(T[i] < sim->T_min){
			if(sim->verbose) cout << "WARNING: resetting T["<<i<<"] from "<<T[i]<<" to "<<sim->T_min<<endl;
			T[i] = sim->T_min;
		}
		if(Ye[i] > sim->Ye_max){
			if(sim->verbose) cout << "WARNING: resetting Ye["<<i<<"] from "<<Ye[i]<<" to "<<sim->Ye_max<<endl;
			Ye[i] = sim->Ye_max;
		}
		if(Ye[i] < sim->Ye_min){
			if(sim->verbose) cout << "WARNING: resetting Ye["<<i<<"] from "<<Ye[i]<<" to "<<sim->Ye_min<<endl;
			Ye[i] = sim->Ye_min;
		}
	}

	// read some parameters
	do_annihilation = lua->scalar<int>("do_annihilation");

	// complain if the grid is obviously not right
	if(rho.size()==0){
		cout << "Error: there are no grid zones." << endl;
		exit(5);
	}

	// calculate integrated quantities to check
	double total_rest_mass = 0.0;
	double total_KE        = 0.0;
	double total_TE        = 0.0;
	double total_hvis      = 0.0;
	double Tbar            = 0.0;
	double Yebar           = 0.0;
	int do_visc = lua->scalar<int>("do_visc");
    #pragma omp parallel for reduction(+:total_rest_mass,total_KE,total_TE,total_hvis,Tbar,Yebar)
	for(size_t z_ind=0;z_ind<rho.size();z_ind++){
		// calculate cell rest mass
		double rest_mass   = zone_rest_mass(z_ind);

		// sanity checks
		PRINT_ASSERT(rho[z_ind],>=,0.0);
		PRINT_ASSERT(zone_com_3volume(z_ind),>=,0.0);
		PRINT_ASSERT(rest_mass,>=,0);

		// calculating totals and averages
		total_rest_mass += rest_mass;
		Tbar            += T[z_ind] * rest_mass;
		Yebar           += Ye[z_ind] * rest_mass;
		total_KE        += (zone_lorentz_factor(z_ind) - 1.0) * rest_mass * pc::c*pc::c;
		total_TE        += rest_mass   / pc::m_n * pc::k * T[z_ind];
		if(do_visc) total_hvis += H_vis[z_ind] * rho[z_ind] * zone_com_3volume(z_ind);
	}

	// write out useful info about the grid
	if (rank0){
		cout << "#   mass = " << total_rest_mass/pc::M_sun << " Msun" << endl;
		cout << "#   <T> = " << Tbar/total_rest_mass*pc::k_MeV << " MeV" << endl;
		cout << "#   <Ye> = " << Yebar/total_rest_mass << endl;
		cout << "#   KE = " << total_KE << " erg" << endl;
		cout << "#   TE = " << total_TE << " erg" << endl;
		if(do_visc) cout << "#   hvis(nonrel) = " << total_hvis << " erg/s" << endl;
	}

    // get the frequency grid
	string neutrino_type = lua->scalar<string>("neutrino_type");
	if(neutrino_type=="NuLib") nulib_get_nu_grid(nu_grid_axis);
	else if(neutrino_type=="GR1D") Neutrino_GR1D::set_nu_grid(lua,&nu_grid_axis);
	else{
        double minval = 0;
        double trash, tmp=0;
        vector<double> bintops = vector<double>(0);
    	int n_nu   = lua->scalar<int>("nugrid_n");
    	if(n_nu>0){
    		double nu_start = lua->scalar<double>("nugrid_start");
    		double nu_stop  = lua->scalar<double>("nugrid_stop");
    		PRINT_ASSERT(nu_stop,>,nu_start);
    		nu_grid_axis = Axis(nu_start/pc::h_MeV,nu_stop/pc::h_MeV,n_nu);
    	}
    	else{
    		string nugrid_filename = lua->scalar<string>("nugrid_filename");
    		ifstream nugrid_file;
    		nugrid_file.open(nugrid_filename.c_str());
    		nugrid_file >> trash >> minval;
    		minval /= pc::h_MeV;
    		vector<double> binmid;
    		while(nugrid_file >> trash >> tmp){
    			tmp /= pc::h_MeV;
    			double last = bintops.size()>0 ? bintops[bintops.size()-1] : minval;
    			PRINT_ASSERT(tmp,>,last);
    			bintops.push_back(tmp);
    			binmid.push_back(0.5 * (last + tmp));
    		}
    		nugrid_file.close();
            nu_grid_axis = Axis(minval,bintops,binmid);
    	}
	}
	PRINT_ASSERT(nu_grid_axis.size(),>,0);

    //==================================//
    // set up the spectrum in each zone //
    //==================================//
    if(rank0) cout << "#   Setting up the distribution function..." << flush;
    string distribution_type = lua->scalar<string>("distribution_type");

    double minval = 0;
    double trash, tmp=0;
    vector<double> bintops = vector<double>(0);
    distribution.resize(insim->species_list.size());

    for(size_t s=0; s<insim->species_list.size(); s++){

    	//-- POLAR SPECTRUM -----------------
    	if(distribution_type == "Polar"){
    		Axis tmp_mugrid, tmp_phigrid;
     		int n_mu = lua->scalar<int>("distribution_nmu");
    		int n_phi = lua->scalar<int>("distribution_nphi");
    		if(n_mu>0 && n_phi>0){
    			tmp_mugrid = Axis(-1,1,n_mu);
    			tmp_phigrid = Axis(-pc::pi, pc::pi, n_phi);
    		}
    		else{
    			// mu ---------------
    			string mugrid_filename = lua->scalar<string>("distribution_mugrid_filename");
    			ifstream mugrid_file;
    			mugrid_file.open(mugrid_filename.c_str());
    			mugrid_file >> trash >> minval;
    			bintops = vector<double>(0);
    			vector<double> binmid;
    			while(mugrid_file >> trash >> tmp){
    				if(bintops.size()>0) PRINT_ASSERT(tmp,>,bintops[bintops.size()-1]);
    				else PRINT_ASSERT(tmp,>,minval);
    				bintops.push_back(tmp);
    				double last = bintops.size()==1 ? minval : bintops[bintops.size()-2];
    				binmid.push_back(0.5 * (last + bintops[bintops.size()-1]));
    			}
    			bintops[bintops.size()-1] = 1.0;
    			PRINT_ASSERT(minval,==,-1);
    			tmp_mugrid = Axis(minval, bintops, binmid);

    			// phi --------------
    			string phigrid_filename = lua->scalar<string>("distribution_phigrid_filename");
    			ifstream phigrid_file;
    			phigrid_file.open(phigrid_filename.c_str());
    			phigrid_file >> trash >> minval;
    			minval -= pc::pi;
    			bintops = vector<double>(0);
    			binmid = vector<double>(0);
    			while(phigrid_file >> trash >> tmp){
    				tmp -= pc::pi;
    				if(bintops.size()>0) PRINT_ASSERT(tmp,>,bintops[bintops.size()-1]);
    				else PRINT_ASSERT(tmp,>,minval);
    				bintops.push_back(tmp);
    				double last = bintops.size()==1 ? minval : bintops[bintops.size()-2];
    				binmid.push_back(0.5 * (last + bintops[bintops.size()-1]));
    			}
    			bintops[bintops.size()-1] = pc::pi;
    			PRINT_ASSERT(minval,==,-pc::pi);
    			tmp_phigrid = Axis(minval, bintops, binmid);
    		}
    		distribution[s] = new PolarSpectrumArray<NDIMS>;
    		((PolarSpectrumArray<NDIMS>*)distribution[s])->init(xAxes,nu_grid_axis, tmp_mugrid, tmp_phigrid);
    	}

    	//-- MOMENT SPECTRUM --------------------
    	else if(distribution_type == "Moments"){
    		distribution[s] = new MomentSpectrumArray<NDIMS>;
    		((MomentSpectrumArray<NDIMS>*)distribution[s])->init(xAxes,nu_grid_axis);
    	}

    	//-- RADIAL MOMENT SPECTRUM --------------------
    	else if(distribution_type == "RadialMoments"){
    		distribution[s] = new RadialMomentSpectrumArray<NDIMS>;
    		((RadialMomentSpectrumArray<NDIMS>*)distribution[s])->init(xAxes,nu_grid_axis);
    	}

    	//-- RADIAL MOMENT SPECTRUM --------------------
    	else if(distribution_type == "GR1D"){
    		distribution[s] = new GR1DSpectrumArray;
    		((GR1DSpectrumArray*)distribution[s])->init(xAxes,nu_grid_axis);
    	}

    	//-- CATCH ------------------
    	else{
    		cout << "Distribution type not found." << endl;
    		assert(0);
    	}

    }

    if(rank0) cout << "finished." << endl << flush;

	// set up the data structures
	abs_opac.resize(sim->species_list.size());
	scat_opac.resize(sim->species_list.size());
	scattering_delta.resize(sim->species_list.size());
	partial_scat_opac.resize(sim->species_list.size());
	spectrum.resize(sim->species_list.size());
	vector<Axis> axes = xAxes;
	if(do_annihilation) fourforce_annihil.set_axes(axes);
	fourforce_abs.set_axes(axes);
	fourforce_emit.set_axes(axes);
	l_abs.set_axes(axes);
	l_emit.set_axes(axes);

	axes.push_back(nu_grid_axis);
	for(size_t s=0; s<sim->species_list.size(); s++){
		abs_opac[s].set_axes(axes);
		scat_opac[s].set_axes(axes);

	    //===========================//
		// intialize output spectrum // only if child didn't
		//===========================//
		Axis tmp_mugrid, tmp_phigrid;
		if(spectrum[s].size()==0){
			int nmu  = lua->scalar<int>("spec_n_mu");
			int nphi = lua->scalar<int>("spec_n_phi");
			tmp_mugrid = Axis(-1,1,nmu);
			tmp_phigrid = Axis(-pc::pi,pc::pi, nphi);
			vector<Axis> dummy_spatial_axes(0);
			spectrum[s].init(dummy_spatial_axes, nu_grid_axis, tmp_mugrid, tmp_phigrid);
		}
		PRINT_ASSERT(spectrum[s].size(),>,0);
	}

	if(sim->use_scattering_kernels==1){
		for(size_t s=0; s<sim->species_list.size(); s++){
			partial_scat_opac[s].resize(nu_grid_axis.size());
			for(size_t igout=0; igout<nu_grid_axis.size(); igout++)
				partial_scat_opac[s][igout].set_axes(axes);
		}
		axes.push_back(nu_grid_axis);
		for(size_t s=0; s<sim->species_list.size(); s++)
			scattering_delta[s].set_axes(axes);
	}
}

void Grid::write_zones(const int iw)
{
	PRINT_ASSERT(rho.size(),>,0);

	// output all zone data in hdf5 format
	string filename = Transport::filename("fluid",iw,".h5");
	H5::H5File file(filename, H5F_ACC_TRUNC);
	file.createGroup("/axes");

	// axes
	for(int dir=0; dir<NDIMS; dir++){
		string datasetname = string("/axes/x")+to_string(dir)+string("(cm)");
		xAxes[dir].write_HDF5(datasetname, file);
	}
	distribution[0]->write_hdf5_coordinates(file, "/axes/distribution");
	nu_grid_axis.write_HDF5("/axes/frequency(Hz)",file);
	spectrum[0].write_hdf5_coordinates(file,"/axes/spectrum");

	// write fluid quantities
	rho.write_HDF5(file,"rho(g|ccm,tet)");
	T.write_HDF5(file,"T_gas(K,tet)");
	Ye.write_HDF5(file,"Ye");
	H_vis.write_HDF5(file,"H_vis(erg|s|g,tet)");
	fourforce_abs.write_HDF5(file,"four-force[abs](erg|ccm|s,tet)");
	fourforce_emit.write_HDF5(file,"four-force[emit](erg|ccm|s,tet)");
	l_abs.write_HDF5(file,"l_abs(1|s|ccm,tet)");
	l_emit.write_HDF5(file,"l_emit(1|s|ccm,tet)");
	if(DO_GR){
		lapse.write_HDF5(file,"lapse");
	}
	if(do_annihilation>0) fourforce_annihil.write_HDF5(file,"annihilation_4force(erg|ccm|s,tet)");
	for(size_t s=0; s<distribution.size(); s++){
		distribution[s]->write_hdf5_data(file, "distribution"+to_string(s)+"(erg|ccm,tet)");
		spectrum[s].write_hdf5_data(file,"spectrum"+to_string(s)+"(erg|s)");
	}

	write_child_zones(file);
}

void Grid::read_zones(string filename)
{
	H5::H5File file(filename, H5F_ACC_RDONLY);

	// axes
	for(int dir=0; dir<NDIMS; dir++){
		string datasetname = string("/axes/x")+to_string(dir)+string("(cm)");
		xAxes[dir].read_HDF5(datasetname, file);
	}
	distribution[0]->read_hdf5_coordinates(file, "/axes/distribution");
	nu_grid_axis.read_HDF5("/axes/frequency(Hz)",file);
	spectrum[0].read_hdf5_coordinates(file,"/axes/spectrum");

	// fluid quantities
	rho.read_HDF5(file,"rho(g|ccm,tet)");
	T.read_HDF5(file,"T_gas(K,tet)");
	Ye.read_HDF5(file,"Ye");
	H_vis.read_HDF5(file,"H_vis(erg|s|g,tet)");
	fourforce_abs.read_HDF5(file,"four-force[abs](erg|ccm|s,tet)");
	fourforce_emit.read_HDF5(file,"four-force[emit](erg|ccm|s,tet)");
	l_abs.read_HDF5(file,"l_abs(1|s|ccm,tet)");
	l_emit.read_HDF5(file,"l_emit(1|s|ccm,tet)");
	if(DO_GR){
		lapse.read_HDF5(file,"lapse");
	}
	if(do_annihilation>0) fourforce_annihil.read_HDF5(file,"annihilation_4force(erg|ccm|s,tet)");
	for(size_t s=0; s<distribution.size(); s++){
		distribution[s]->read_hdf5_data(file, "distribution"+to_string(s)+"(erg|ccm,tet)");
		spectrum[s].read_hdf5_data(file,"spectrum"+to_string(s)+"(erg|s)");
	}

	read_child_zones(file);
}

double Grid::zone_rest_mass(const int z_ind) const{
	return rho[z_ind] * zone_com_3volume(z_ind);
}

double Grid::zone_com_3volume(const int z_ind) const{
	// assumes v is orthonormal in cm/s
	PRINT_ASSERT(z_ind,>=,0);
	double result = zone_lab_3volume(z_ind) * zone_lorentz_factor(z_ind);
	PRINT_ASSERT(result,>,0);
	return result;
}


double Grid::total_rest_mass() const{
	double mass = 0;
	#pragma omp parallel for reduction(+:mass)
	for(size_t z_ind=0; z_ind<rho.size(); z_ind++) mass += zone_rest_mass(z_ind);
	return mass;
}

// radius given coordinates
double Grid::zone_4volume(const int z_ind) const{
	return zone_lab_3volume(z_ind) * (DO_GR ? lapse[z_ind] : 1.0);
}

void Grid::interpolate_metric(EinsteinHelper *eh) const{
  assert(DO_GR);

  // first, the lapse
  eh->g.alpha = lapse.interpolate(eh->icube_vol);
  PRINT_ASSERT(eh->g.alpha,>,0);

  // second, the shift and three-metric
  eh->g.betaup = interpolate_shift(*eh);
  eh->g.gammalow.data = interpolate_3metric(*eh);

  // fill in the rest of the metric values
  eh->g.update();
}
