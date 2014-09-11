#include <mpi.h>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
#include "grid_general.h"
#include "Lua.h"
#include "global_options.h"

//------------------------------------------------------------
// initialize the grid
//------------------------------------------------------------
void grid_general::init(Lua* lua)
{
	// read the model file or fill in custom model
	std::string model_file = lua->scalar<std::string>("model_file");
	if(model_file == "custom") custom_model(lua);
	else read_model_file(lua);

	// complain if the grid is obviously not right
	if(z.size()==0){
		cout << "Error: there are no grid zones." << endl;
		exit(5);
	}
}


//------------------------------------------------------------
// Write the grid information out to a file
//------------------------------------------------------------
void grid_general::write_zones(const int iw) const
{
	assert(z.size()>0);


	ofstream outf;
	transport::open_file("fluid",iw,outf);
	write_header(outf);
	vector<int> dir_ind;
	for (unsigned z_ind=0; z_ind<z.size(); z_ind++)
	{
		zone_directional_indices(z_ind, dir_ind);
		if(dir_ind.size()>0) if(dir_ind[dir_ind.size()-1]==0) outf << endl;
		write_line(outf,z_ind);
	}
	outf.close();
}

double grid_general::zone_rest_mass(const int z_ind) const{
	return z[z_ind].rho*zone_lab_volume(z_ind)*transport::lorentz_factor(z[z_ind].v);
}

double grid_general::zone_comoving_volume(const int z_ind) const{
	// assumes v is orthonormal in cm/s
	vector<double> v_tmp(3,0);
	for(unsigned i=0; i<z[z_ind].v.size(); i++) v_tmp[i] = z[z_ind].v[i];
	if(!good_zone(z_ind)) return 0;
	else return zone_lab_volume(z_ind) * transport::lorentz_factor(v_tmp);
}

void grid_general::write_header(ofstream& outf) const{
	outf << setprecision(4);
	outf << scientific;
	outf << "# ";
	vector<double> r;
	zone_coordinates(0,r);
	const int dimensionality = r.size();
	for(int i=0; i<dimensionality; i++) outf << "r[" << i << "] ";
	outf << "e_rad(erg/ccm)  rho(g/ccm)  T_gas(MeV)  Ye  t_therm  t_lep  |v|(cm/s)  H-C(erg/g/s)  dYe_dt(1/s)" << endl;
}

void grid_general::write_line(ofstream& outf, const int z_ind) const{
	vector<double> r;
	zone_coordinates(z_ind,r);
	assert(r.size()>=0);

	for(unsigned i=0; i<r.size(); i++) outf << r[i] << " ";

	outf << z[z_ind].e_rad << " ";
	outf << z[z_ind].rho   << " ";
	outf << z[z_ind].T*pc::k_MeV << " ";
	outf << z[z_ind].Ye    << " ";

	outf << 1.0 / fabs(1.0/z[z_ind].t_eabs - 1.0/z[z_ind].t_eemit) << " ";
	outf << 1.0 / fabs(1.0/z[z_ind].t_labs - 1.0/z[z_ind].t_lemit) << " ";

	outf << zone_speed2(z_ind) << " ";

	double net_neutrino_energy_source = (z[z_ind].e_abs - z[z_ind].e_emit) / z[z_ind].rho - z[z_ind].H_com;
	outf << net_neutrino_energy_source << " ";

	double n_baryons_per_ccm = z[z_ind].rho / transport::mean_mass(z[z_ind].Ye);
	double dYe_dt = (z[z_ind].l_abs - z[z_ind].l_emit) / n_baryons_per_ccm;
	outf << dYe_dt << " ";

	outf << endl;
}

bool grid_general::good_zone(const int z_ind) const{
	//const zone* z = &(grid->z[z_ind]);
	//return (z->rho >= rho_min && z->rho <= rho_max &&
	//  	    z->Ye  >= Ye_min  && z->Ye  <=  Ye_max &&
	//        z->T   >=  T_min  && z->T   <=   T_max);
	return zone_speed2(z_ind) < pc::c*pc::c;
}



//------------------------------------
// get the velocity squared of a zone
//------------------------------------
double grid_general::zone_speed2(const int z_ind) const{
	assert(z_ind >= 0);
	assert(z_ind < (int)z.size());
	double speed2 = transport::dot(z[z_ind].v,z[z_ind].v);
	return speed2;
}
