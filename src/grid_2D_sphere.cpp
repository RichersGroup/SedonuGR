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
#include <iomanip>
#include "global_options.h"
#include "grid_2D_sphere.h"
#include "transport.h"

using namespace std;
namespace pc = physical_constants;

//------------------------------------------------------------
// initialize the zone geometry from model file
//------------------------------------------------------------
void grid_2D_sphere::read_model_file(Lua* lua)
{
	std::string model_type = lua->scalar<std::string>("model_type");
	if(model_type == "flash") read_flash_file(lua);
	else if(model_type == "custom") custom_model(lua);
	else{
		cout << "ERROR: model type unknown." << endl;
		exit(8);
	}
}

void grid_2D_sphere::read_flash_file(Lua* lua)
{
	// verbocity
	int my_rank;
	MPI_Comm_rank( MPI_COMM_WORLD, &my_rank );
	const int rank0 = (my_rank == 0);

	// generalHDF5 variables
	H5::DataSet dataset;
	H5::DataSpace space;
	H5::CompType comptype;

	// open the model files
	if(rank0) cout << "# Reading the model files..." << endl;
	string model_filename   = lua->scalar<string>("model_file"  );
	string xCoords_filename = lua->scalar<string>("xCoords_file");
	string yCoords_filename = lua->scalar<string>("yCoords_file");
	H5::H5File file(model_filename, H5F_ACC_RDONLY);
	ifstream xCoords_file, yCoords_file;
	xCoords_file.open(xCoords_filename.c_str());
	yCoords_file.open(yCoords_filename.c_str());
	if(xCoords_file.fail()){
		if(rank0) cout << "Error: can't read the xCoords file." << xCoords_filename << endl;
		exit(4);
	}
	if(yCoords_file.fail()){
		if(rank0) cout << "Error: can't read the xCoords file." << yCoords_filename << endl;
		exit(4);
	}

	//====================================================//
	// get the general properties from "/integer scalars" //
	//====================================================//

	// get the database, etc
	dataset = file.openDataSet("/integer scalars");
	space = dataset.getSpace();
	comptype = dataset.getCompType();

	// set up the structures in memory
	const size_t stringsize = 81;
	struct pair_t{
		char name[stringsize];
		int value;
	};
	H5::StrType string_type(H5T_STRING,stringsize);
	H5::CompType mempair_t(sizeof(pair_t));
	mempair_t.insertMember("name",HOFFSET(pair_t,name),string_type);
	mempair_t.insertMember("value",HOFFSET(pair_t,value),H5::PredType::STD_I32LE);

	// make sure we have the right dataset
	hsize_t dim;
	const int list_rank = 1;
	PRINT_ASSERT(space.getSimpleExtentNdims(),==,list_rank);                 // 1D array
	space.getSimpleExtentDims(&dim);
	PRINT_ASSERT(dataset.getTypeClass(),==,H5T_COMPOUND);                  // filled with structs
	PRINT_ASSERT(comptype.getNmembers(),==,2);                               // each struct has 2 elements
	PRINT_ASSERT(comptype.getMemberName(0),==,"name");                       // first element is the name
	PRINT_ASSERT(comptype.getMemberName(1),==,"value");                      // second element is the value
	PRINT_ASSERT(comptype.getMemberClass(0),==,H5T_STRING);                // the type of the first element is a string
	PRINT_ASSERT(comptype.getMemberClass(1),==,H5T_INTEGER);               // the type of the second element is a 32bit little endian int
	PRINT_ASSERT(comptype.getMemberDataType(0).getSize(),==,stringsize-1); // the name is 80 characters long, 1 extra for the null terminate

	// read the data
	vector<pair_t> integer_data(dim);
	dataset.read(&(integer_data[0]),mempair_t);
	for(unsigned i=0; i<dim; i++) integer_data[i].name[stringsize-1] = '\0';
	PRINT_ASSERT(trim(string(integer_data[0].name)),==,string("nxb")); // make sure we're looking at the right fields
	PRINT_ASSERT(trim(string(integer_data[1].name)),==,string("nyb"));
	PRINT_ASSERT(trim(string(integer_data[2].name)),==,string("nzb"));
	PRINT_ASSERT(trim(string(integer_data[3].name)),==,string("dimensionality"));
	PRINT_ASSERT(trim(string(integer_data[4].name)),==,string("iprocs"));
	PRINT_ASSERT(trim(string(integer_data[5].name)),==,string("jprocs"));
	PRINT_ASSERT(trim(string(integer_data[6].name)),==,string("kprocs"));
	const int nxb = integer_data[0].value;
	const int nyb = integer_data[1].value;
	PRINT_ASSERT(integer_data[2].value,==,1); // 2D dataset should have 1 thickness in the z direction
	PRINT_ASSERT(integer_data[3].value,==,2); // 2D dataset
	const int iprocs = integer_data[4].value;
	const int jprocs = integer_data[5].value;
	PRINT_ASSERT(integer_data[6].value,==,1); // 2D dataset cannot be split in z direction

	// deduce the global structure
	const int nr     = nxb*iprocs;
	const int ntheta = nyb*jprocs;
	const int n_zones = nr*ntheta;
	z.resize(n_zones);


	//=========================//
	// read in the actual data //
	//=========================//

	// check that everything makes sense with one of the datasets
	const int dataset_rank = 4;
	hsize_t dims[dataset_rank];
	dataset = file.openDataSet("/dens");
	space = dataset.getSpace();
	space.getSimpleExtentDims(dims);
	PRINT_ASSERT(dataset.getTypeClass(),==,H5T_FLOAT);
	PRINT_ASSERT(space.getSimpleExtentNdims(),==,dataset_rank);
	PRINT_ASSERT((int)dims[0],==,iprocs*jprocs);
	PRINT_ASSERT((int)dims[1],==,1);
	PRINT_ASSERT((int)dims[2],==,nyb);
	PRINT_ASSERT((int)dims[3],==,nxb);

	// read the data
	float dens[dims[0]][dims[1]][dims[2]][dims[3]]; // g/ccm
	float velx[dims[0]][dims[1]][dims[2]][dims[3]]; // cm/s
	float vely[dims[0]][dims[1]][dims[2]][dims[3]]; // cm/s
	float angz[dims[0]][dims[1]][dims[2]][dims[3]]; // cm^2/s
	float efrc[dims[0]][dims[1]][dims[2]][dims[3]]; //
	float temp[dims[0]][dims[1]][dims[2]][dims[3]]; // K
	float hvis[dims[0]][dims[1]][dims[2]][dims[3]]; // erg/g/s (viscous heating)
	float gamn[dims[0]][dims[1]][dims[2]][dims[3]]; // 1/s     (net rate of change of Ye)
	float ncfn[dims[0]][dims[1]][dims[2]][dims[3]]; // erg/g/s (net     charged-current neutrino heating-cooling)
	float nprs[dims[0]][dims[1]][dims[2]][dims[3]]; // erg/g/s (net non-charged-current neutrino heating-cooling)
	float eint[dims[0]][dims[1]][dims[2]][dims[3]]; // erg/g   (internal energy density)
	float atms[dims[0]][dims[1]][dims[2]][dims[3]]; //         (hydrogen mass fraction)
	float neut[dims[0]][dims[1]][dims[2]][dims[3]]; //         (neutron  mass fraction)
	float prot[dims[0]][dims[1]][dims[2]][dims[3]]; //         (proton   mass fraction)
	float alfa[dims[0]][dims[1]][dims[2]][dims[3]]; //         (alpha    mass fraction)
	dataset = file.openDataSet("/dens");
	dataset.read(&(dens[0][0][0][0]),H5::PredType::IEEE_F32LE);
	dataset = file.openDataSet("/velx");
	dataset.read(&(velx[0][0][0][0]),H5::PredType::IEEE_F32LE);
	dataset = file.openDataSet("/vely");
	dataset.read(&(vely[0][0][0][0]),H5::PredType::IEEE_F32LE);
	dataset = file.openDataSet("/angz");
	dataset.read(&(angz[0][0][0][0]),H5::PredType::IEEE_F32LE);
	dataset = file.openDataSet("/efrc");
	dataset.read(&(efrc[0][0][0][0]),H5::PredType::IEEE_F32LE);
	dataset = file.openDataSet("/temp");
	dataset.read(&(temp[0][0][0][0]),H5::PredType::IEEE_F32LE);
	dataset = file.openDataSet("/hvis");
	dataset.read(&(hvis[0][0][0][0]),H5::PredType::IEEE_F32LE);
	dataset = file.openDataSet("/gamn");
	dataset.read(&(gamn[0][0][0][0]),H5::PredType::IEEE_F32LE);
	dataset = file.openDataSet("/ncfn");
	dataset.read(&(ncfn[0][0][0][0]),H5::PredType::IEEE_F32LE);
	dataset = file.openDataSet("/nprs");
	dataset.read(&(nprs[0][0][0][0]),H5::PredType::IEEE_F32LE);
	dataset = file.openDataSet("/eint");
	dataset.read(&(eint[0][0][0][0]),H5::PredType::IEEE_F32LE);
	dataset = file.openDataSet("/atms");
	dataset.read(&(atms[0][0][0][0]),H5::PredType::IEEE_F32LE);
	dataset = file.openDataSet("/neut");
	dataset.read(&(neut[0][0][0][0]),H5::PredType::IEEE_F32LE);
	dataset = file.openDataSet("/prot");
	dataset.read(&(prot[0][0][0][0]),H5::PredType::IEEE_F32LE);
	dataset = file.openDataSet("/alfa");
	dataset.read(&(alfa[0][0][0][0]),H5::PredType::IEEE_F32LE);
	dataset.close();
	file.close();


	//=========================//
	// read in the coordinates //
	//=========================//

	// check that the files have the correct number of lines
	const int nghost = 4;
	int x_linecount=0, y_linecount=0;
	string line;
	float trash = 0;
	while(getline(xCoords_file,line)) x_linecount++;
	while(getline(yCoords_file,line)) y_linecount++;
	PRINT_ASSERT(x_linecount,==,iprocs*nxb+2*nghost);
	PRINT_ASSERT(y_linecount,==,jprocs*nyb+2*nghost);

	// read x (r) coordinates
	r_out.resize(nr);
	xCoords_file.clear();                                   // clear bad state from EOF
	xCoords_file.seekg(0);                                  // seek back to beginning of the file
	for(int i=0; i<nghost; i++) getline(xCoords_file,line); // trash the first four lines (ghost points)
	for(int i=0; i<nr; i++){
		if(i==0) xCoords_file >> r_out.min;
		else xCoords_file >> trash;
		xCoords_file >> trash;
		xCoords_file >> r_out[i];
		xCoords_file >> trash;
	}

	// read the y (theta) coordinates
	theta_out.resize(ntheta);
	yCoords_file.clear();  // clear bad state from EOF
	yCoords_file.seekg(0); // seek back to beginning of the file
	for(int i=0; i<nghost; i++) getline(yCoords_file,line); // trash the first four lines (ghost points)
	for(int i=0; i<ntheta; i++){
		if(i==0) yCoords_file >> theta_out.min;
		else yCoords_file >> trash;
		yCoords_file >> trash;
		yCoords_file >> theta_out[i];
		yCoords_file >> trash;
	}


	//===============//
	// fill the grid //
	//===============//
	int do_visc = lua->scalar<int>("do_visc");
	const int kb = 0;
	double newtonian_eint_total = 0;
	double net_HmC = 0;
	double net_dyedt = 0;
	double newtonian_mass = 0;
	double Mnu = 0;
	const double gamma_max = 2.0;
	const double speed_max = pc::c * sqrt(1.0 - 1.0/gamma_max);
    #pragma omp parallel for collapse(3) reduction(+:newtonian_eint_total,net_HmC,net_dyedt,newtonian_mass,Mnu)
	for(unsigned proc=0; proc<dims[0]; proc++)
		for(unsigned jb=0; jb<dims[2]; jb++)
			for(unsigned ib=0; ib<dims[3]; ib++){
				// indices. moving by one proc in the x direction increases proc by 1
				const int i_global = (proc%iprocs)*nxb + ib;
				const int j_global = (proc/iprocs)*nyb + jb;
				const int z_ind = zone_index(i_global, j_global);
				PRINT_ASSERT(i_global,<,nr);
				PRINT_ASSERT(j_global,<,ntheta);
				PRINT_ASSERT(z_ind,<,n_zones);

				// zone position
				double r[2];
				zone_coordinates(z_ind,r,2);

				// zone values
				z[z_ind].rho               = dens[proc][kb][jb][ib];
				z[z_ind].T                 = temp[proc][kb][jb][ib];
				z[z_ind].Ye                = efrc[proc][kb][jb][ib];
				if(do_visc) z[z_ind].H_vis = hvis[proc][kb][jb][ib];
				double vr              = velx[proc][kb][jb][ib];
				double vtheta          = vely[proc][kb][jb][ib];
				double vphi            = angz[proc][kb][jb][ib]/r[0]/sin(r[1]);
				double speed = sqrt(vr*vr + vtheta*vtheta + vphi*vphi);
				if(speed > speed_max){
					vr     *= speed_max / speed;
					vtheta *= speed_max / speed;
					vphi   *= speed_max / speed;
					if(rank0) cout << "WARNING: velocity of superluminal cell at {r,theta}={" << r[0] << "," << r[1] << "} set to gamma=" << gamma_max << endl;
				}
				PRINT_ASSERT(ARRSIZE(z[z_ind].u),==,3);
				z[z_ind].u[0] = vr;
				z[z_ind].u[1] = vtheta;
				z[z_ind].u[2] = vphi;
				PRINT_ASSERT(zone_speed2(z_ind),<,pc::c*pc::c);
				PRINT_ASSERT(z[z_ind].rho,>=,0.0);
				PRINT_ASSERT(z[z_ind].T,>=,0.0);
				PRINT_ASSERT(z[z_ind].Ye,>=,0.0);
				PRINT_ASSERT(z[z_ind].Ye,<=,1.0);
				newtonian_eint_total += eint[proc][kb][jb][ib]             * z[z_ind].rho * zone_lab_volume(z_ind);
				net_HmC += (ncfn[proc][kb][jb][ib]-nprs[proc][kb][jb][ib]) * z[z_ind].rho * zone_lab_volume(z_ind);
				net_dyedt += gamn[proc][kb][jb][ib] * z[z_ind].rho * zone_lab_volume(z_ind);
				newtonian_mass += z[z_ind].rho * zone_lab_volume(z_ind);
				if(z[z_ind].H_vis < ncfn[proc][kb][jb][ib]) Mnu += z[z_ind].rho * zone_lab_volume(z_ind);
	}
	if(rank0){
		cout << "#   Newtonian total internal energy: " << newtonian_eint_total << " erg" << endl;
		cout << "#   H-C (Newtonian, neutrino only): " << net_HmC << " erg/s" << endl;
		cout << "#   dYe_dt (Newtonian): " << net_dyedt/newtonian_mass << " 1/s" << endl;
		cout << "#   M_nu (Newtonian): " << Mnu/physical_constants::M_sun << " Msun" << endl;
	}


	//================================//
	// output timescales in data file //
	//================================//
	if(rank0){
		int width = 15;

		vector<double> z_gamn(z.size());
		vector<double> z_ncfn(z.size());
		vector<double> z_nprs(z.size());
		vector<double> z_eint(z.size());
		for(unsigned proc=0; proc<dims[0]; proc++)
			for(unsigned jb=0; jb<dims[2]; jb++)
				for(unsigned ib=0; ib<dims[3]; ib++){
					// indices. moving by one proc in the x direction increases proc by 1
					const int i_global = (proc%iprocs)*nxb + ib;
					const int j_global = (proc/iprocs)*nyb + jb;
					const int z_ind = zone_index(i_global, j_global);
					PRINT_ASSERT(i_global,<,nr);
					PRINT_ASSERT(j_global,<,ntheta);
					PRINT_ASSERT(z_ind,<,n_zones);

					z_gamn[z_ind] = gamn[proc][kb][jb][ib];
					z_ncfn[z_ind] = ncfn[proc][kb][jb][ib];
					z_nprs[z_ind] = nprs[proc][kb][jb][ib];
					z_eint[z_ind] = eint[proc][kb][jb][ib];
				}

		ofstream outf;
		outf.open("rates_flash.dat");
		outf << setw(width) << "r(cm)";
		outf << setw(width) << "theta";
		outf << setw(width) << "gamn(1/s)";
		outf << setw(width) << "ncfn(erg/g/s)";
		outf << setw(width) << "nprs(erg/g/s)";
		outf << setw(width) << "hvis(erg/g/s)";
		outf << setw(width) << "eint(erg/g)";
		outf << setw(width) << endl;
		for(unsigned z_ind=0; z_ind<z.size(); z_ind++){
			// zone position
			double r[2];
			zone_coordinates(z_ind,r,2);

			//double gamma = transport::lorentz_factor(z[z_ind].v);
			//double m_zone = z[z_ind].rho * zone_lab_volume(z_ind)*gamma;
			//double t_lep = 1.0/z_gamn[z_ind];
			//double t_therm = m_zone*pc::k*z[z_ind].T/z_ncfn[z_ind];
			if(z_ind%theta_out.size() == 0) outf << endl;
			outf << setw(width) << r[0];
			outf << setw(width) << r[1];
			outf << setw(width) << z_gamn[z_ind];
			outf << setw(width) << z_ncfn[z_ind];
			outf << setw(width) << z_nprs[z_ind];
			outf << setw(width) << z[z_ind].H_vis;
			outf << setw(width) << z_eint[z_ind];
			outf << endl;
		}
		outf.close();

		//=================================//
		// output composition in data file //
		//=================================//
		vector<double> z_atms(z.size());
		vector<double> z_neut(z.size());
		vector<double> z_prot(z.size());
		vector<double> z_alfa(z.size());
		for(unsigned proc=0; proc<dims[0]; proc++)
			for(unsigned jb=0; jb<dims[2]; jb++)
				for(unsigned ib=0; ib<dims[3]; ib++){
					// indices. moving by one proc in the x direction increases proc by 1
					const int i_global = (proc%iprocs)*nxb + ib;
					const int j_global = (proc/iprocs)*nyb + jb;
					const int z_ind = zone_index(i_global, j_global);
					PRINT_ASSERT(i_global,<,nr);
					PRINT_ASSERT(j_global,<,ntheta);
					PRINT_ASSERT(z_ind,<,n_zones);

					z_atms[z_ind] = atms[proc][kb][jb][ib];
					z_neut[z_ind] = neut[proc][kb][jb][ib];
					z_prot[z_ind] = prot[proc][kb][jb][ib];
					z_alfa[z_ind] = alfa[proc][kb][jb][ib];
				}

		outf.open("initial_composition.dat");
		outf << setw(width) << "r(cm)";
		outf << setw(width) << "theta";
		outf << setw(width) << "x_hydrogen";
		outf << setw(width) << "x_neutrons";
		outf << setw(width) << "x_protons";
		outf << setw(width) << "x_alpha";
		outf << setw(width) << "(1-x_total)";
		outf << setw(width) << "derived_Ye";
		outf << setw(width) << "data_Ye";
		outf << setw(width) << "T(MeV)";
		outf << endl;
		for(unsigned z_ind=0; z_ind<z.size(); z_ind++){
			// zone position
			double r[2];
			zone_coordinates(z_ind,r,2);

			// calculate electron fraction
			double abar = 1.0 / (z_prot[z_ind] + z_neut[z_ind] + z_alfa[z_ind]/4.0);
			double zbar = abar * (z_prot[z_ind] + 2.0*z_alfa[z_ind]/4.0);
			double ye_calculated  = zbar/abar;//total_protons / total_nucleons;

			if(z_ind%theta_out.size() == 0) outf << endl;
			outf << setw(width) << r[0];
			outf << setw(width) << r[1];
			outf << setw(width) << z_atms[z_ind];
			outf << setw(width) << z_neut[z_ind];
			outf << setw(width) << z_prot[z_ind];
			outf << setw(width) << z_alfa[z_ind];
			outf << setw(width) << 1.0-(z_atms[z_ind]+z_neut[z_ind]+z_prot[z_ind]+z_alfa[z_ind]);
			outf << setw(width) << ye_calculated;
			outf << setw(width) << z[z_ind].Ye;
			outf << setw(width) << z[z_ind].T*pc::k_MeV;
			outf << endl;
		}
		outf.close();
	}
}

//------------------------------------------------------------
// Write a custom model here if you like
//------------------------------------------------------------
void grid_2D_sphere::custom_model(Lua* lua)
{
	// verbocity
	int my_rank;
	MPI_Comm_rank( MPI_COMM_WORLD, &my_rank );
	const int rank0 = (my_rank == 0);
	if(rank0) cout << "#   Reading 1D model file, mapping to 2D" << endl;

	// open up the model file, complaining if it fails to open
	string model_file = "neutron_star.mod";
	ifstream infile;
	infile.open(model_file.c_str());
	if(infile.fail()){
		cout << "Error: can't read the model file." << model_file << endl;
		exit(4);
	}

	// geometry of model
	infile >> grid_type;
	if(grid_type != "1D_sphere"){
		cout << "Error: grid_type parameter disagrees with the model file." << endl;
	}

	// number of zones
	int r_zones,theta_zones;
	theta_zones = 4;
	infile >> r_zones;
	PRINT_ASSERT(r_zones,>,0);
	z.resize(r_zones*theta_zones);
	r_out.resize(r_zones);
	theta_out.resize(theta_zones);

	// read zone properties
	infile >> r_out.min;
	PRINT_ASSERT(r_out.min,>=,0);
	theta_out.min = 0;
	for(int j=0; j<theta_zones; j++) theta_out[j] = theta_out.min + (double)(j+1) * pc::pi / (double)theta_zones;
	for(int i=0; i<r_zones; i++)
	{
		infile >> r_out[i];
		PRINT_ASSERT(r_out[i],>,(i==0 ? r_out.min : r_out[i-1]));

		int base_ind = zone_index(i,0);
		infile >> z[base_ind].rho;
		infile >> z[base_ind].T;
		infile >> z[base_ind].Ye;
		z[base_ind].H_vis = 0;
		PRINT_ASSERT(ARRSIZE(z[base_ind].u),==,3);
		z[base_ind].u[0] = 0;
		z[base_ind].u[1] = 0;
		z[base_ind].u[2] = 0;
		PRINT_ASSERT(z[base_ind].rho,>=,0);
		PRINT_ASSERT(z[base_ind].T,>=,0);
		PRINT_ASSERT(z[base_ind].Ye,>=,0);
		PRINT_ASSERT(z[base_ind].Ye,<=,1.0);

		for(int j=0; j<theta_zones; j++){
			int z_ind = zone_index(i,j);
			z[z_ind] = z[base_ind];
		}
	}

	infile.close();
}

//------------------------------------------------------------
// Return the zone index containing the position x
//------------------------------------------------------------
int grid_2D_sphere::zone_index(const double x[3], const int xsize) const
{
	PRINT_ASSERT(xsize,==,3);
	double r  = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
	double theta = atan2(sqrt(x[0]*x[0] + x[1]*x[1]), x[2]);
	PRINT_ASSERT(r,>=,0);
	PRINT_ASSERT(theta,>=,0);
	if(fabs(theta-pc::pi)<tiny) theta = pc::pi-tiny;
	PRINT_ASSERT(theta,<=,pc::pi);

	// check if off the boundaries
	if(r     <  r_out.min                    ) return -1;
	if(r     >= r_out[r_out.size()-1]        ) return -2;
	if(theta <  theta_out.min                ) return -2;
	if(theta >= theta_out[theta_out.size()-1]) return -2;

	// find in zone array using stl algorithm upper_bound and subtracting iterators
	const int i =     r_out.locate(r    );
	const int j = theta_out.locate(theta);
	const int z_ind = zone_index(i,j);
	PRINT_ASSERT(z_ind,>=,0);
	PRINT_ASSERT(z_ind,<,(int)z.size());
	return z_ind;
}

//----------------------------------------------------------------
// Return the zone index corresponding to the directional indices
//----------------------------------------------------------------
int grid_2D_sphere::zone_index(const int i, const int j) const
{
	PRINT_ASSERT(i,>=,0);
	PRINT_ASSERT(j,>=,0);
	PRINT_ASSERT(i,<,(int)r_out.size());
	PRINT_ASSERT(j,<,(int)theta_out.size());
	const int z_ind = i*theta_out.size() + j;
	PRINT_ASSERT(z_ind,<,(int)z.size());
	return z_ind;
}


//------------------------------------------------------------
// return volume of zone
//------------------------------------------------------------
double grid_2D_sphere::zone_lab_volume(const int z_ind) const
{
	PRINT_ASSERT(z_ind,>=,0);
	PRINT_ASSERT(z_ind,<,(int)z.size());
	int dir_ind[2];
	zone_directional_indices(z_ind,dir_ind,2);
	const unsigned i = dir_ind[0];
	const unsigned j = dir_ind[1];
	const double r0     =     r_out.bottom(i);
	const double theta0 = theta_out.bottom(j);
	const double r1     =     r_out[i];
	const double theta1 = theta_out[j];
	const double vol = 2.0*pc::pi/3.0 * (cos(theta0) - cos(theta1)) * (r1*r1*r1 - r0*r0*r0);
	PRINT_ASSERT(vol,>=,0);
	return vol;
}


//------------------------------------------------------------
// return length of zone
//------------------------------------------------------------
double grid_2D_sphere::zone_min_length(const int z_ind) const
{
	int dir_ind[2];
	zone_directional_indices(z_ind,dir_ind,2);
	const unsigned i = dir_ind[0];
	const unsigned j = dir_ind[1];

	// the 'minimum lengts' are just approximate.
	const double r_len     = (    r_out[i] -     r_out.bottom(i));
	const double theta_len = (theta_out[j] - theta_out.bottom(j)) * r_out.bottom(i);

	// if r_in is zero, there will be problems, but simulations would not have done this.
	if(r_out.bottom(i) == 0) return r_len;
	else return min(r_len, theta_len);

}

//------------------------------------------------------------
// Return the cell-center spherical coordinates of the cell
//------------------------------------------------------------
void grid_2D_sphere::zone_coordinates(const int z_ind, double r[2], const int rsize) const
{
	PRINT_ASSERT(z_ind,>=,0);
	PRINT_ASSERT(z_ind,<,(int)(r_out.size()*theta_out.size()));
	PRINT_ASSERT(rsize,==,2);

	int dir_ind[2];
	zone_directional_indices(z_ind, dir_ind, 2);
	const unsigned i = dir_ind[0];
	const unsigned j = dir_ind[1];

	const double r0     =     r_out.bottom(i);
	const double theta0 = theta_out.bottom(j);
	r[0] = 0.5 * (r0     +     r_out[i]);
	r[1] = 0.5 * (theta0 + theta_out[j]);
}


//-------------------------------------------
// get directional indices from zone index
//-------------------------------------------
void grid_2D_sphere::zone_directional_indices(const int z_ind, int dir_ind[2], const int size) const
{
	PRINT_ASSERT(z_ind,>=,0);
	PRINT_ASSERT(z_ind,<,(int)z.size());
	PRINT_ASSERT(size,==,2);
	dir_ind[0] = z_ind / theta_out.size(); // r index
	dir_ind[1] = z_ind % theta_out.size(); // theta index
	PRINT_ASSERT(dir_ind[0],>=,0);
	PRINT_ASSERT(dir_ind[1],>=,0);
	PRINT_ASSERT(dir_ind[0],<,(int)    r_out.size());
	PRINT_ASSERT(dir_ind[1],<,(int)theta_out.size());
}


//------------------------------------------------------------
// sample a random cartesian position within the spherical shell
//------------------------------------------------------------
void grid_2D_sphere::cartesian_sample_in_zone(const int z_ind, const double rand[3], const int randsize, double x[3], const int xsize) const
{
	PRINT_ASSERT(z_ind,>=,0);
	PRINT_ASSERT(z_ind,<,(int)z.size());
	PRINT_ASSERT(randsize,==,3);
	PRINT_ASSERT(xsize,==,3);

	// radius and theta indices
	int dir_ind[2];
	zone_directional_indices(z_ind,dir_ind,2);
	int i = dir_ind[0];
	int j = dir_ind[1];

	// inner and outer coordinates of shell
	double  r0 =         r_out.bottom(i);
	double mu0 = cos(theta_out.bottom(j));
	double  r1 =         r_out[i];
	double mu1 = cos(theta_out[j]);

	// sample radial position in shell using a probability integral transform
	double radius = pow( rand[0]*(r1*r1*r1 - r0*r0*r0) + r0*r0*r0, 1./3.);
	PRINT_ASSERT(radius,>=,r0*(1.-tiny));
	PRINT_ASSERT(radius,<=,r1*(1.+tiny));
	radius = max(r0,radius);
	radius = min(r1,radius);

	// sample cos(theta) uniformily
	double mu = mu0 + (mu1-mu0)*rand[1];
	mu = max(mu0, mu);
	mu = min(mu1, mu);
	double sin_theta = sqrt(1-mu*mu);
	PRINT_ASSERT(sin_theta,>=,-1.0);
	PRINT_ASSERT(sin_theta,<=,1.0);

	// sample phi uniformily
	double phi = 2.0*pc::pi*rand[2];

	// set the real 3-d coordinates. remember, z is along the symmetry axis
	x[0] = radius*sin_theta*cos(phi);
	x[1] = radius*sin_theta*sin(phi);
	x[2] = radius*mu;
}


//------------------------------------------------------------
// get the cartesian velocity vector (cm/s)
//------------------------------------------------------------
void grid_2D_sphere::cartesian_velocity_vector(const double x[3], const int xsize, double v[3], const int vsize, int z_ind) const
{
	PRINT_ASSERT(xsize,==,3);
	PRINT_ASSERT(vsize,==,3);
	if(z_ind < 0) z_ind = zone_index(x,xsize);
	PRINT_ASSERT(z_ind,>=,-1);

	// if within inner sphere, z_ind=-1. Leave velocity at 0.
	if(z_ind >= 0){

		// radius in zone
		double r    = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
		double rhat = sqrt(x[0]*x[0] + x[1]*x[1]);
		int along_axis = (rhat/r < tiny);

		// Based on position, calculate what the 3-velocity is
		PRINT_ASSERT(ARRSIZE(z[z_ind].u),==,3);
		double vr     = z[z_ind].u[0];
		double vtheta = z[z_ind].u[1];
		double vphi   = z[z_ind].u[2];

		double vr_cart[3];
		vr_cart[0] = vr * x[0]/r;
		vr_cart[1] = vr * x[1]/r;
		vr_cart[2] = vr * x[2]/r;

		double vtheta_cart[3];
		vtheta_cart[0] =  (along_axis ? 0 : vtheta * x[2]/r * x[0]/rhat );
		vtheta_cart[1] =  (along_axis ? 0 : vtheta * x[2]/r * x[1]/rhat );
		vtheta_cart[2] = -vtheta * rhat/r;

		double vphi_cart[3];
		vphi_cart[0] = (along_axis ? 0 : -vphi * x[1]/rhat );
		vphi_cart[1] = (along_axis ? 0 :  vphi * x[0]/rhat );
		vphi_cart[2] = 0;

		// remember, symmetry axis is along the z-axis
		for(int i=0; i<3; i++) v[i] = vr_cart[i] + vtheta_cart[i] + vphi_cart[i];

		// check for pathological case
		if (r == 0){ // set everything to 0
			v[0] = 0;
			v[1] = 0;
			v[2] = 0;
		}
	}
}


//------------------------------------------------------------
// Write the grid information out to a file
//------------------------------------------------------------
void grid_2D_sphere::write_rays(int iw) const
{
	PRINT_ASSERT(iw,>=,0);
	ofstream outf;
	unsigned i=0,j=0;
	double r[2];
	string filename = "";

	// along theta=0
	filename = transport::filename("ray_t0",iw,".dat");
	outf.open(filename.c_str());
	write_header(outf);
	j = 0;
	for(i=0; i<r_out.size(); i++){
		int z_ind = zone_index(i,j);
		zone_coordinates(z_ind,r,2);
		write_line(outf,z_ind);
	}
	outf.close();

	// along theta=pi/2
	filename = transport::filename("ray_t.5",iw,".dat");
	outf.open(filename.c_str());
	write_header(outf);
	j = theta_out.size()/2;
	for(i=0; i<r_out.size(); i++){
		int z_ind = zone_index(i,j);
		zone_coordinates(z_ind,r,2);
		write_line(outf,z_ind);
	}
	outf.close();

	// along theta=pi
	filename = transport::filename("ray_t1",iw,".dat");
	outf.open(filename.c_str());
	write_header(outf);
	j = theta_out.size()-1;
	for(i=0; i<r_out.size(); i++){
		int z_ind = zone_index(i,j);
		zone_coordinates(z_ind,r,2);
		write_line(outf,z_ind);
	}
	outf.close();

	// along theta
	filename = transport::filename("ray_r.5",iw,".dat");
	outf.open(filename.c_str());
	write_header(outf);
	i = r_out.size()/2;
	for(j=0; j<theta_out.size(); j++){
		int z_ind = zone_index(i,j);
		zone_coordinates(z_ind,r,2);
		write_line(outf,z_ind);
	}
	outf.close();
}


//------------------------------------------------------------
// Reflect off the outer boundary
//------------------------------------------------------------
void grid_2D_sphere::reflect_outer(particle *p) const{
	PRINT_ASSERT(r_out.size(),>=,1);
	double r0 = (r_out.size()==1 ? r_out.min : r_out.size()-2);
	double dr = r_out[r_out.size()-1] - r0;
	PRINT_ASSERT( fabs(p->r() - r_out[r_out.size()-1]),<,tiny*dr);
	double velDotRhat = p->mu();

	// invert the radial component of the velocity
	for(int i=0; i<3; i++) p->D[i] -= 2.*velDotRhat * p->x[i]/p->r();
	transport::normalize(p->D,3);

	// put the particle just inside the boundary
	double newR = r_out[r_out.size()-1] - tiny*dr;
	for(int i=0; i<3; i++) p->x[i] = p->x[i]/p->r()*newR;

	// must be inside the boundary, or will get flagged as escaped
	PRINT_ASSERT(zone_index(p->x,3),>=,0);
}

//------------------------------------------------------------
// Reflect off the symmetry boundaries
//------------------------------------------------------------
void grid_2D_sphere::symmetry_boundaries(particle *p) const{
// not implemented - does nothing
}

//------------------------------------------------------------
// Find distance to outer boundary
//------------------------------------------------------------
double grid_2D_sphere::lab_dist_to_boundary(const particle *p) const{
	// Theta = angle between radius vector and direction (Pi if outgoing)
	// Phi   = Pi - Theta (angle on the triangle) (0 if outgoing)
	double Rout  = r_out[r_out.size()-1];
	double Rin   = r_out.min;
	double r  = p->r();
	double mu = p->mu();
	double d_outer_boundary = numeric_limits<double>::infinity();
	double d_inner_boundary = numeric_limits<double>::infinity();
	PRINT_ASSERT(r,<,Rout);
	PRINT_ASSERT(zone_index(p->x,3),>=,-1);

	// distance to inner boundary
	if(r >= Rin){
		double radical = r*r*(mu*mu-1.0) + Rin*Rin;
		if(Rin>0 && mu<0 && radical>=0){
			d_inner_boundary = -r*mu - sqrt(radical);
			PRINT_ASSERT(d_inner_boundary,<=,sqrt(Rout*Rout-Rin*Rin)*(1.0+tiny));
		}
	}
	else{
		d_inner_boundary = -r*mu + sqrt(r*r*(mu*mu-1.0) + Rin*Rin);
		PRINT_ASSERT(d_inner_boundary,<=,2.*Rin);
	}
	if(d_inner_boundary<0 && fabs(d_inner_boundary/Rin)<tiny*(r_out[0]-Rin)) d_inner_boundary = tiny*(r_out[0]-Rin);
	PRINT_ASSERT(d_inner_boundary,>=,0);

	// distance to outer boundary
	d_outer_boundary = -r*mu + sqrt(r*r*(mu*mu-1.0) + Rout*Rout);
	if(d_outer_boundary<=0 && fabs(d_outer_boundary/Rin)<tiny*(Rout-r_out[r_out.size()-1])) d_outer_boundary = tiny*(Rout-r_out[r_out.size()-1]);
	PRINT_ASSERT(d_outer_boundary,>,0);
	PRINT_ASSERT(d_outer_boundary,<=,2.*Rout);

	// distances to the theta boundaries - NOT IMPLEMENTED THETA BOUNDARIES
	PRINT_ASSERT( fabs((theta_out[theta_out.size()-1] - theta_out.min) - pc::pi),<,tiny);
	double theta_dist = INFINITY;

	// make sure the particle ends up in a reasonable place
	const double r_dist = min(d_inner_boundary, d_outer_boundary);
	return min(r_dist,theta_dist);
}


double grid_2D_sphere::zone_radius(const int z_ind) const{
	PRINT_ASSERT(z_ind,>=,0);
	PRINT_ASSERT(z_ind,<,(int)z.size());

	// radius and theta indices
	int dir_ind[2];
	zone_directional_indices(z_ind,dir_ind,2);
	int i = dir_ind[0];

	return r_out[i];
}

//-----------------------------
// Dimensions of the grid
//-----------------------------
void grid_2D_sphere::dims(hsize_t dims[2], const int size) const{
	PRINT_ASSERT(size,==,dimensionality());
	dims[0] = r_out.size();
	dims[1] = theta_out.size();
}

//----------------------------------------------------
// Write the coordinates of the grid points to the hdf5 file
//----------------------------------------------------
void grid_2D_sphere::write_hdf5_coordinates(H5::H5File file) const
{
	// useful quantities
	H5::DataSet dataset;
	H5::DataSpace dataspace;

	// get dimensions
	hsize_t coord_dims[2];
	dims(coord_dims,2);
	for(unsigned i=0; i<2; i++) coord_dims[i]++; //make room for min value

	// write r coordinates
	float tmp0[coord_dims[0]];
	dataspace = H5::DataSpace(1,&coord_dims[0]);
	dataset = file.createDataSet("grid_r(cm)",H5::PredType::IEEE_F32LE,dataspace);
	tmp0[0] = r_out.min;
	for(unsigned i=1; i<r_out.size()+1; i++) tmp0[i] = r_out[i-1];
	dataset.write(&tmp0[0],H5::PredType::IEEE_F32LE);
	dataset.close();

	// write theta coordinates
	float tmp1[coord_dims[1]];
	dataspace = H5::DataSpace(1,&coord_dims[1]);
	dataset = file.createDataSet("grid_theta(radians)",H5::PredType::IEEE_F32LE,dataspace);
	tmp1[0] = theta_out.min;
	for(unsigned i=1; i<theta_out.size()+1; i++) tmp1[i] = theta_out[i-1];
	dataset.write(&tmp1[0],H5::PredType::IEEE_F32LE);
	dataset.close();
}
