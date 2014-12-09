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
	// MPI stuff
	int my_rank,n_procs;
	MPI_Comm_rank( MPI_COMM_WORLD, &my_rank );
	MPI_Comm_size( MPI_COMM_WORLD, &n_procs);
	bool rank0 = (my_rank==0);

	// read the model file or fill in custom model
	std::string model_file = lua->scalar<std::string>("model_file");
	int do_relativity = lua->scalar<int>("do_relativity");
	if(model_file == "custom") custom_model(lua);
	else read_model_file(lua);

	// complain if the grid is obviously not right
	if(z.size()==0){
		cout << "Error: there are no grid zones." << endl;
		exit(5);
	}

	// calculate integrated quantities to check
	double total_nonrel_mass = 0.0;
	double total_rest_mass   = 0.0;
	double total_rel_KE      = 0.0;
	double total_nonrel_KE   = 0.0;
	double total_rel_TE      = 0.0;
	double total_nonrel_TE   = 0.0;
	double total_hvis        = 0.0;
	int do_visc = lua->scalar<int>("do_visc");
    #pragma omp parallel for reduction(+:total_nonrel_mass, total_rest_mass, total_rel_KE, total_nonrel_KE, total_rel_TE, total_nonrel_TE, total_hvis)
	for(unsigned z_ind=0;z_ind<z.size();z_ind++){
		double rest_mass   = z[z_ind].rho * zone_comoving_volume(z_ind);
		assert(rest_mass >= 0);
		double nonrel_mass = z[z_ind].rho * zone_lab_volume(z_ind);
		assert(nonrel_mass >= 0);
		vector<double> r;
		zone_coordinates(z_ind,r);

		//if(grid->z[z_ind].rho > 1.0e8){ // && r[1] > pc::pi/3.0 && r[1] < pc::pi/2.0){
		total_rest_mass += rest_mass;
		total_nonrel_mass += nonrel_mass;
		total_rel_KE    += (rest_mass>0 ? (transport::lorentz_factor(z[z_ind].v) - 1.0) * rest_mass * pc::c*pc::c : 0);
		total_nonrel_KE += 0.5 * nonrel_mass * zone_speed2(z_ind);
		total_rel_TE    += (rest_mass>0 ? rest_mass   / pc::m_n * pc::k * z[z_ind].T : 0);
		total_nonrel_TE += nonrel_mass / pc::m_n * pc::k * z[z_ind].T;
		if(do_visc) total_hvis += z[z_ind].H_vis * z[z_ind].rho * zone_comoving_volume(z_ind);
		if(!do_relativity) for(unsigned i=0; i<z[z_ind].v.size(); i++) z[z_ind].v[i] = 0;
		//}
	}
	if (rank0){
		cout << "#   mass = " << total_rest_mass << " g (nonrel: " << total_nonrel_mass << " g)" <<endl;
		cout << "#   KE = " << total_rel_KE << " erg (nonrel: " << total_nonrel_KE << " erg)" << endl;
		cout << "#   TE = " << total_rel_TE << " erg (nonrel: " << total_nonrel_TE << " erg)" << endl;
		if(do_visc) cout << "#   hvis = " << total_hvis << " erg/s" << endl;
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
	unsigned c=0;
	for(unsigned i=0; i<r.size(); i++) outf << ++c << "-r[" << i << "]  ";
	outf << ++c << "-comoving_volume(ccm)  ";
	outf << ++c << "-rho(g/ccm)  ";
	outf << ++c << "-T_gas(MeV)  ";
	outf << ++c << "-Ye  ";
	outf << ++c << "-e_rad(erg/ccm)  ";
	outf << ++c << "-H_vis(erg/s/g)  ";
	//outf << ++c << "-t_therm  ";
	//outf << ++c << "-t_lep  ";
	outf << ++c << "-|v|(cm/s)  ";
	outf << ++c << "-H-C(erg/g/s)  ";
	outf << ++c << "-dYe_dt(1/s)  ";
	outf << ++c << "-annihilation_rate(erg/ccm/s)  ";
	if(z[0].distribution.size()>0) outf << ++c << "-integrated_edens(erg/ccm)  ";
	for(unsigned s=0; s<z[0].distribution.size(); s++){
		for(unsigned g=0; g<z[0].distribution[s].size(); g++){
			outf << ++c << "-s"<<s<<"g"<<g<<"edens(erg/ccm)  ";
		}
	}
	outf << endl;
}

void grid_general::write_line(ofstream& outf, const int z_ind) const{
	vector<double> r;
	zone_coordinates(z_ind,r);

	for(unsigned i=0; i<r.size(); i++) outf << r[i] << " ";

	outf << zone_comoving_volume(z_ind) << "\t";
	outf << z[z_ind].rho   << "\t";
	outf << z[z_ind].T*pc::k_MeV << "\t";
	outf << z[z_ind].Ye    << "\t";
	outf << z[z_ind].e_rad << "\t";
	outf << z[z_ind].H_vis << "\t";

	//outf << 1.0 / fabs(1.0/z[z_ind].t_eabs - 1.0/z[z_ind].t_eemit) << "\t";
	//outf << 1.0 / fabs(1.0/z[z_ind].t_labs - 1.0/z[z_ind].t_lemit) << "\t";

	outf << zone_speed2(z_ind) << "\t";

	double net_neutrino_energy_source = (z[z_ind].e_abs - z[z_ind].e_emit) / z[z_ind].rho - z[z_ind].H_vis;
	outf << net_neutrino_energy_source << "\t";

	double n_baryons_per_ccm = z[z_ind].rho / transport::mean_mass(z[z_ind].Ye);
	double dYe_dt = (z[z_ind].nue_abs-z[z_ind].anue_abs - z[z_ind].l_emit) / n_baryons_per_ccm;
	outf << dYe_dt << "\t";

	outf << z[z_ind].Q_annihil << "\t";

	double integrated_edens = 0;
	for(unsigned s=0; s<z[z_ind].distribution.size(); s++){
		vector<double> direction_integrated_edens;
		z[z_ind].distribution[s].integrate_over_direction(direction_integrated_edens);
		assert(direction_integrated_edens.size() == z[z_ind].distribution[s].size());
		for(unsigned g=0; g<z[z_ind].distribution[s].size(); g++){
			outf << direction_integrated_edens[g] << "\t";
			integrated_edens += direction_integrated_edens[g];
		}
	}
	if(z[z_ind].distribution.size()>0) outf << integrated_edens << "\t";

	outf << endl;
}

bool grid_general::good_zone(const int z_ind) const{
	//const zone* z = &(grid->z[z_ind]);
	//return (z->rho >= rho_min && z->rho <= rho_max &&
	//  	    z->Ye  >= Ye_min  && z->Ye  <=  Ye_max &&
	//        z->T   >=  T_min  && z->T   <=   T_max);
	if(z_ind < 0) return false;
	else if(zone_speed2(z_ind) > pc::c*pc::c) return false;
	else return true;
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
