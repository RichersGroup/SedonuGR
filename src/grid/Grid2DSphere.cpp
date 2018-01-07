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
#include "Grid2DSphere.h"
#include "Transport.h"
#include "H5Cpp.h"

using namespace std;
namespace pc = physical_constants;

Grid2DSphere::Grid2DSphere(){
	PRINT_ASSERT(NDIMS,==,2);
	grid_type = "Grid2DSphere";
	tetrad_rotation = spherical;
}

//------------------------------------------------------------
// initialize the zone geometry from model file
//------------------------------------------------------------
void Grid2DSphere::read_model_file(Lua* lua)
{
	std::string model_type = lua->scalar<std::string>("model_type");
	if(model_type == "Flash") read_flash_file(lua);
	else if(model_type == "Nagakura") read_nagakura_file(lua);
	else{
		cout << "ERROR: model type unknown." << endl;
		exit(8);
	}

	vr.calculate_slopes();
	vtheta.calculate_slopes();
	vphi.calculate_slopes();
}
void Grid2DSphere::read_nagakura_file(Lua* lua)
{
	// verbocity
	int MPI_myID;
	MPI_Comm_rank( MPI_COMM_WORLD, &MPI_myID );
	const int rank0 = (MPI_myID == 0);
	double trash=0, minval=0, tmp=0;
	vector<double> bintops, binmid;

	// open the model files
	if(rank0) cout << "# Reading the model file..." << endl;

	//=========================//
	// read in the radial grid //
	//=========================//
	string rgrid_filename = lua->scalar<string>("Grid2DSphere_Nagakura_rgrid_file");
	ifstream rgrid_file;
	rgrid_file.open(rgrid_filename.c_str());
	rgrid_file >> trash >> minval;
	bintops = vector<double>(0);
	while(rgrid_file >> trash >> tmp){
		if(bintops.size()>0) PRINT_ASSERT(tmp,>,bintops[bintops.size()-1]);
		else PRINT_ASSERT(tmp,>,minval);
		bintops.push_back(tmp);
		double last = bintops.size()==1 ? minval : bintops[bintops.size()-2];
		binmid.push_back(0.5 * (last + bintops[bintops.size()-1]));
	}
	rAxis = Axis(minval, bintops, binmid);
	int nr = rAxis.size();

	//========================//
	// read in the theta grid //
	//========================//
	string thetagrid_filename = lua->scalar<string>("Grid2DSphere_Nagakura_thetagrid_file");
	ifstream thetagrid_file;
	vector<double> tmp_vector(1,0);
	thetagrid_file.open(thetagrid_filename.c_str());
	thetagrid_file >> trash >> tmp_vector[0];
	while(thetagrid_file >> trash >> tmp) tmp_vector.push_back(tmp);

	// reverse the vector
	int ntheta = tmp_vector.size()-1;
	bintops = vector<double>(ntheta,0);
	binmid = vector<double>(ntheta,0);
	minval = 0;
	for(int i=0; i<ntheta; i++){
		bintops[i] = tmp_vector[ntheta-1-i];
		double last = i==0 ? minval : bintops[i-1];
		binmid[i] = 0.5 * (last + bintops[i]);
		PRINT_ASSERT(bintops[i],>,last);
	}
	PRINT_ASSERT(fabs(bintops[ntheta-1]-pc::pi),<,TINY);
	bintops[ntheta-1] = pc::pi;
	binmid[ntheta-1] = 0.5 * (bintops[ntheta-1] + bintops[ntheta-2]);
	thetaAxis = Axis(minval, bintops, binmid);

	vector<Axis> axes;
	axis_vector(axes);
	rho.set_axes(axes);
	T.set_axes(axes);
	Ye.set_axes(axes);
	H_vis.set_axes(axes);
	vr.set_axes(axes);
	vtheta.set_axes(axes);
	vphi.set_axes(axes);

	// write grid properties
	if(rank0) cout << "#   nr=" << nr << "\trmin=" << rAxis.min << "\trmax=" << rAxis.top[nr-1] << endl;
	if(rank0) cout << "#   nt=" << ntheta << "\ttmin=" << thetaAxis.min << "\ttmax=" << thetaAxis.top[ntheta-1] << endl;

	//===========================//
	// read the fluid properties //
	//===========================//
	string model_file = lua->scalar<string>("model_file");
	ifstream infile;
	infile.open(model_file.c_str());
	if(infile.fail()){
		if(rank0) cout << "Error: can't read the model file." << model_file << endl;
		exit(4);
	}
	for(int itheta=ntheta-1; itheta>=0; itheta--) // Hiroki's theta is backwards
		for(int ir=0; ir<nr; ir++){
			unsigned z_ind = zone_index(ir,itheta);
			double trash;

			// read the contents of a single line
			infile >> trash; // r sin(theta)
			infile >> trash; // r cos(theta)
			infile >> rho[z_ind]; // g/ccm
			infile >> Ye[z_ind];
			infile >> T[z_ind]; // MeV
			infile >> vr[z_ind]; // cm/s
			infile >> vtheta[z_ind]; // 1/s

			// skip the rest of the line
			for(int k=9; k<=165; k++) infile >> trash;

			// convert units
			//z[z_ind].u[1] *= r_out.center(ir);
			//z[z_ind].u[2] *= r_out.center(ir);
			T[z_ind] /= pc::k_MeV;

			// sanity checks
			PRINT_ASSERT(rho[z_ind],>=,0.0);
			PRINT_ASSERT(T[z_ind],>=,0.0);
			PRINT_ASSERT(Ye[z_ind],>=,0.0);
			PRINT_ASSERT(Ye[z_ind],<=,1.0);
		}
}

void Grid2DSphere::read_flash_file(Lua* lua)
{
	// verbocity
	int MPI_myID;
	MPI_Comm_rank( MPI_COMM_WORLD, &MPI_myID );
	const int rank0 = (MPI_myID == 0);

	// generalHDF5 variables
	H5::DataSet dataset;
	H5::DataSpace space;
	H5::CompType comptype;

	// open the model files
	if(rank0) cout << "# Reading the model files..." << endl;
	string model_filename   = lua->scalar<string>("model_file"  );
	string xCoords_filename = lua->scalar<string>("Grid2DSphere_Flash_xCoords_file");
	string yCoords_filename = lua->scalar<string>("Grid2DSphere_Flash_yCoords_file");
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
	vector<double> bintop(nr), binmid(nr);
	double minval;
	xCoords_file.clear();                                   // clear bad state from EOF
	xCoords_file.seekg(0);                                  // seek back to beginning of the file
	for(int i=0; i<nghost; i++) getline(xCoords_file,line); // trash the first four lines (ghost points)
	for(int i=0; i<nr; i++){
		if(i==0) xCoords_file >> minval;
		else xCoords_file >> trash;
		xCoords_file >> trash;
		xCoords_file >> bintop[i];
		xCoords_file >> trash;
		double last = i==0 ? minval : bintop[i-1];
		binmid[i] = 0.5 * (last + bintop[i]);
	}
	rAxis = Axis(minval, bintop, binmid);

	// read the y (theta) coordinates
	bintop.resize(ntheta);
	binmid.resize(ntheta);
	yCoords_file.clear();  // clear bad state from EOF
	yCoords_file.seekg(0); // seek back to beginning of the file
	for(int i=0; i<nghost; i++) getline(yCoords_file,line); // trash the first four lines (ghost points)
	for(int i=0; i<ntheta; i++){
		if(i==0) yCoords_file >> minval;
		else yCoords_file >> trash;
		yCoords_file >> trash;
		yCoords_file >> bintop[i];
		yCoords_file >> trash;
		double last = i==0 ? minval : bintop[i-1];
		binmid[i] = 0.5 * (last + bintop[i]);
	}
	thetaAxis = Axis(minval, bintop, binmid);

	vector<Axis> axes;
	axis_vector(axes);
	rho.set_axes(axes);
	T.set_axes(axes);
	Ye.set_axes(axes);
	H_vis.set_axes(axes);
	vr.set_axes(axes);
	vtheta.set_axes(axes);
	vphi.set_axes(axes);

	//===============//
	// fill the grid //
	//===============//
	int do_visc = lua->scalar<int>("do_visc");
	const int kb = 0;
	const double gamma_max = 2.0;
	const double speed_max = pc::c * sqrt(1.0 - 1.0/gamma_max);
    #pragma omp parallel for collapse(3)
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
				rho[z_ind]               = dens[proc][kb][jb][ib];
				T[z_ind]                 = temp[proc][kb][jb][ib];
				Ye[z_ind]                = efrc[proc][kb][jb][ib];
				if(do_visc) H_vis[z_ind] = hvis[proc][kb][jb][ib];
				double tmp_vr            = velx[proc][kb][jb][ib];
				double tmp_vtheta        = vely[proc][kb][jb][ib];
				double tmp_vphi          = angz[proc][kb][jb][ib]/r[0]/sin(r[1]);
				double speed = sqrt(tmp_vr*tmp_vr + tmp_vtheta*tmp_vtheta + tmp_vphi*tmp_vphi);
				if(speed > speed_max){
					tmp_vr     *= speed_max / speed;
					tmp_vtheta *= speed_max / speed;
					tmp_vphi   *= speed_max / speed;
					if(rank0) cout << "WARNING: velocity of superluminal cell at {r,theta}={" << r[0] << "," << r[1] << "} set to gamma=" << gamma_max << endl;
				}
				vr[z_ind] = tmp_vr;
				vtheta[z_ind] = tmp_vtheta;
				vphi[z_ind] = tmp_vphi;
				PRINT_ASSERT(rho[z_ind],>=,0.0);
				PRINT_ASSERT(T[z_ind],>=,0.0);
				PRINT_ASSERT(Ye[z_ind],>=,0.0);
				PRINT_ASSERT(Ye[z_ind],<=,1.0);
	}


	//================================//
	// output timescales in data file //
	//================================//
	if(rank0){
		int width = 15;

		vector<double> z_gamn(rho.size());
		vector<double> z_ncfn(rho.size());
		vector<double> z_nprs(rho.size());
		vector<double> z_eint(rho.size());
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
		for(unsigned z_ind=0; z_ind<rho.size(); z_ind++){
			// zone position
			double r[2];
			zone_coordinates(z_ind,r,2);

			//double gamma = transport::lorentz_factor(z[z_ind].v);
			//double m_zone = z[z_ind].rho * zone_lab_volume(z_ind)*gamma;
			//double t_lep = 1.0/z_gamn[z_ind];
			//double t_therm = m_zone*pc::k*z[z_ind].T/z_ncfn[z_ind];
			if(z_ind%thetaAxis.size() == 0) outf << endl;
			outf << setw(width) << r[0];
			outf << setw(width) << r[1];
			outf << setw(width) << z_gamn[z_ind];
			outf << setw(width) << z_ncfn[z_ind];
			outf << setw(width) << z_nprs[z_ind];
			outf << setw(width) << H_vis[z_ind];
			outf << setw(width) << z_eint[z_ind];
			outf << endl;
		}
		outf.close();

		//=================================//
		// output composition in data file //
		//=================================//
		vector<double> z_atms(rho.size());
		vector<double> z_neut(rho.size());
		vector<double> z_prot(rho.size());
		vector<double> z_alfa(rho.size());
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
		for(unsigned z_ind=0; z_ind<rho.size(); z_ind++){
			// zone position
			double r[2];
			zone_coordinates(z_ind,r,2);

			// calculate electron fraction
			double abar = 1.0 / (z_prot[z_ind] + z_neut[z_ind] + z_alfa[z_ind]/4.0);
			double zbar = abar * (z_prot[z_ind] + 2.0*z_alfa[z_ind]/4.0);
			double ye_calculated  = zbar/abar;//total_protons / total_nucleons;

			if(z_ind%thetaAxis.size() == 0) outf << endl;
			outf << setw(width) << r[0];
			outf << setw(width) << r[1];
			outf << setw(width) << z_atms[z_ind];
			outf << setw(width) << z_neut[z_ind];
			outf << setw(width) << z_prot[z_ind];
			outf << setw(width) << z_alfa[z_ind];
			outf << setw(width) << 1.0-(z_atms[z_ind]+z_neut[z_ind]+z_prot[z_ind]+z_alfa[z_ind]);
			outf << setw(width) << ye_calculated;
			outf << setw(width) << Ye[z_ind];
			outf << setw(width) << T[z_ind]*pc::k_MeV;
			outf << endl;
		}
		outf.close();
	}
}

double Grid2DSphere_theta(const double x[3]){
	return atan2(sqrt(x[0]*x[0] + x[1]*x[1]), x[2]);
}

//------------------------------------------------------------
// Return the zone index containing the position x
//------------------------------------------------------------
int Grid2DSphere::zone_index(const double x[3]) const
{
	double r  = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
	double theta = Grid2DSphere_theta(x);
	PRINT_ASSERT(r,>=,0);
	PRINT_ASSERT(theta,>=,0);
	if(fabs(theta-pc::pi)<TINY) theta = pc::pi-TINY;
	PRINT_ASSERT(theta,<=,pc::pi);

	// check if off the boundaries
	if(r     <  rAxis.min                    ) return -1;
	if(r     >= rAxis.top[rAxis.size()-1]        ) return -2;
	if(theta <  thetaAxis.min                ) return -2;
	if(theta >= thetaAxis.top[thetaAxis.size()-1]) return -2;

	// find in zone array using stl algorithm upper_bound and subtracting iterators
	const int i =     rAxis.bin(r    );
	const int j = thetaAxis.bin(theta);
	const int z_ind = zone_index(i,j);
	PRINT_ASSERT(z_ind,>=,0);
	PRINT_ASSERT(z_ind,<,(int)rho.size());
	return z_ind;
}

//----------------------------------------------------------------
// Return the zone index corresponding to the directional indices
//----------------------------------------------------------------
int Grid2DSphere::zone_index(const int i, const int j) const
{
	PRINT_ASSERT(i,>=,0);
	PRINT_ASSERT(j,>=,0);
	PRINT_ASSERT(i,<,(int)rAxis.size());
	PRINT_ASSERT(j,<,(int)thetaAxis.size());
	const int z_ind = i*thetaAxis.size() + j;
	PRINT_ASSERT(z_ind,<,(int)rho.size());
	return z_ind;
}


//------------------------------------------------------------
// return volume of zone
//------------------------------------------------------------
double Grid2DSphere::zone_lab_3volume(const int z_ind) const
{
	PRINT_ASSERT(DO_GR,==,false); // need to include sqrt(detg3)
	PRINT_ASSERT(z_ind,>=,0);
	PRINT_ASSERT(z_ind,<,(int)rho.size());
	vector<unsigned> dir_ind(2);
	zone_directional_indices(z_ind,dir_ind);
	const unsigned i = dir_ind[0];
	const unsigned j = dir_ind[1];
	const double r0     =     rAxis.bottom(i);
	const double theta0 = thetaAxis.bottom(j);
	const double r1     =     rAxis.top[i];
	const double theta1 = thetaAxis.top[j];
	const double vol = 2.0*pc::pi/3.0 * (cos(theta0) - cos(theta1)) * (r1*r1*r1 - r0*r0*r0);
	PRINT_ASSERT(vol,>=,0);
	return vol;
}


// returning 0 causes the min distance to take over in propagate.cpp::which_event

double Grid2DSphere::d_boundary(const EinsteinHelper* eh) const{
	const double r = radius(eh->p.xup);
	PRINT_ASSERT(r,<=,rAxis.top[eh->z_ind]);
	PRINT_ASSERT(r,>=,rAxis.bottom(eh->z_ind));
	const double x = eh->p.xup[0];
	const double y = eh->p.xup[1];
	const double z = eh->p.xup[2];
	const double rp = sqrt(x*x + y*y);
	const double theta = Grid2DSphere_theta(eh->p.xup);

	// get component of k in the radial direction
	double kr = eh->g.dot<4>(eh->e[0],eh->p.kup);
	double dlambda_r = INFINITY;
	if(kr>0) dlambda_r = (rAxis.top[eh->dir_ind[0]] - r   ) / kr;
	if(kr<0) dlambda_r = (r - rAxis.bottom(eh->dir_ind[0])) / kr;
	dlambda_r = abs(dlambda_r);

	// get component of k in the theta direction
	double kt = eh->g.dot<4>(eh->e[1],eh->p.kup);
	double dlambda_t = INFINITY;
	if(kt>0) dlambda_t = r * (thetaAxis.top[eh->dir_ind[1]]  - theta  ) / kt;
	if(kt<0) dlambda_t = r * (theta - thetaAxis.bottom(eh->dir_ind[1])) / kt;
	dlambda_t = abs(dlambda_t);

	double ds_com = eh->kup_tet[3] * min(dlambda_t, dlambda_r);
	return ds_com;
}

//------------------------------------------------------------
// return length of zone
//------------------------------------------------------------
double Grid2DSphere::zone_min_length(const int z_ind) const
{
	vector<unsigned> dir_ind(2);
	zone_directional_indices(z_ind,dir_ind);
	const unsigned i = dir_ind[0];
	const unsigned j = dir_ind[1];

	// the 'minimum lengts' are just approximate.
	const double r_len     = (    rAxis.top[i] -     rAxis.bottom(i));
	const double theta_len = sin(thetaAxis.top[j] - thetaAxis.bottom(j)) * rAxis.bottom(i);

	// if r_in is zero, there will be problems, but simulations would not have done this.
	return min(r_len, theta_len);
}

//------------------------------------------------------------
// Return the cell-center spherical coordinates of the cell
//------------------------------------------------------------
void Grid2DSphere::zone_coordinates(const int z_ind, double r[2], const int rsize) const
{
	PRINT_ASSERT(z_ind,>=,0);
	PRINT_ASSERT(z_ind,<,(int)(rAxis.size()*thetaAxis.size()));
	PRINT_ASSERT(rsize,==,2);

	vector<unsigned> dir_ind(2);
	zone_directional_indices(z_ind, dir_ind);
	const unsigned i = dir_ind[0];
	const unsigned j = dir_ind[1];

	const double r0     =     rAxis.bottom(i);
	const double theta0 = thetaAxis.bottom(j);
	r[0] = 0.5 * (r0     +     rAxis.top[i]);
	r[1] = 0.5 * (theta0 + thetaAxis.top[j]);
}


//-------------------------------------------
// get directional indices from zone index
//-------------------------------------------
void Grid2DSphere::zone_directional_indices(const int z_ind, vector<unsigned>& dir_ind) const
{
	PRINT_ASSERT(z_ind,>=,0);
	PRINT_ASSERT(z_ind,<,(int)rho.size());
	PRINT_ASSERT(dir_ind.size(),==,2);
	dir_ind[0] = z_ind / thetaAxis.size(); // r index
	dir_ind[1] = z_ind % thetaAxis.size(); // theta index
	PRINT_ASSERT(dir_ind[0],>=,0);
	PRINT_ASSERT(dir_ind[1],>=,0);
	PRINT_ASSERT(dir_ind[0],<,(int)    rAxis.size());
	PRINT_ASSERT(dir_ind[1],<,(int)thetaAxis.size());
}


//------------------------------------------------------------
// sample a random cartesian position within the spherical shell
//------------------------------------------------------------
void Grid2DSphere::sample_in_zone(const int z_ind, ThreadRNG* rangen, double x[3]) const
{
	PRINT_ASSERT(z_ind,>=,0);
	PRINT_ASSERT(z_ind,<,(int)rho.size());

	double rand[3];
	rand[0] = rangen->uniform();
	rand[1] = rangen->uniform();
	rand[2] = rangen->uniform();

	// radius and theta indices
	vector<unsigned> dir_ind(2);
	zone_directional_indices(z_ind,dir_ind);
	int i = dir_ind[0];
	int j = dir_ind[1];

	// inner and outer coordinates of shell
	double  r0 =         rAxis.bottom(i);
	double mu0 = cos(thetaAxis.top[j]);
	double  r1 =         rAxis.top[i];
	double mu1 = cos(thetaAxis.bottom(j));
	PRINT_ASSERT(mu1,>,mu0);

	// sample radial position in shell using a probability integral transform
	double radius = pow( rand[0]*(r1*r1*r1 - r0*r0*r0) + r0*r0*r0, 1./3.);
	PRINT_ASSERT(radius,>=,r0*(1.-TINY));
	PRINT_ASSERT(radius,<=,r1*(1.+TINY));
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
void Grid2DSphere::interpolate_fluid_velocity(const double x[3], double v[3], const unsigned dir_ind[NDIMS]) const
{
	// radius in zone
	double r    = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
	double rhat = sqrt(x[0]*x[0] + x[1]*x[1]);
	int along_axis = (rhat/r < TINY);
	double theta = pc::pi/2.0 - atan2(x[2], rhat);
	theta = max(0.0,theta);
	theta = min(pc::pi, theta);

	// Based on position, calculate what the 3-velocity is
	double xvec[2] = {r,theta};
	double tmp_vr     = vr.interpolate(xvec,dir_ind);
	double tmp_vtheta = vtheta.interpolate(xvec,dir_ind);
	double tmp_vphi   = vphi.interpolate(xvec,dir_ind);

	double vr_cart[3];
	vr_cart[0] = tmp_vr * x[0]/r;
	vr_cart[1] = tmp_vr * x[1]/r;
	vr_cart[2] = tmp_vr * x[2]/r;

	double vtheta_cart[3];
	vtheta_cart[0] =  (along_axis ? 0 : tmp_vtheta * x[2]/r * x[0]/rhat );
	vtheta_cart[1] =  (along_axis ? 0 : tmp_vtheta * x[2]/r * x[1]/rhat );
	vtheta_cart[2] = -tmp_vtheta * rhat/r;

	double vphi_cart[3];
	vphi_cart[0] = (along_axis ? 0 : -tmp_vphi * x[1]/rhat );
	vphi_cart[1] = (along_axis ? 0 :  tmp_vphi * x[0]/rhat );
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


//------------------------------------------------------------
// Reflect off the symmetry boundaries
//------------------------------------------------------------
void Grid2DSphere::symmetry_boundaries(EinsteinHelper *lh) const{
// not implemented - does nothing
}

double Grid2DSphere::zone_radius(const int z_ind) const{
	PRINT_ASSERT(z_ind,>=,0);
	PRINT_ASSERT(z_ind,<,(int)rho.size());

	// radius and theta indices
	vector<unsigned> dir_ind(2);
	zone_directional_indices(z_ind,dir_ind);
	int i = dir_ind[0];

	return rAxis.top[i];
}

//-----------------------------
// Dimensions of the grid
//-----------------------------
void Grid2DSphere::dims(hsize_t dims[2], const int size) const{
	PRINT_ASSERT(size,==,(int)dimensionality());
	dims[0] = rAxis.size();
	dims[1] = thetaAxis.size();
}

//----------------------------------------------------
// Write the coordinates of the grid points to the hdf5 file
//----------------------------------------------------
void Grid2DSphere::write_hdf5_coordinates(H5::H5File file) const
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
	tmp0[0] = rAxis.min;
	for(unsigned i=1; i<rAxis.size()+1; i++) tmp0[i] = rAxis.top[i-1];
	dataset.write(&tmp0[0],H5::PredType::IEEE_F32LE);
	dataset.close();

	// write theta coordinates
	float tmp1[coord_dims[1]];
	dataspace = H5::DataSpace(1,&coord_dims[1]);
	dataset = file.createDataSet("grid_theta(radians)",H5::PredType::IEEE_F32LE,dataspace);
	tmp1[0] = thetaAxis.min;
	for(unsigned i=1; i<thetaAxis.size()+1; i++) tmp1[i] = thetaAxis.top[i-1];
	dataset.write(&tmp1[0],H5::PredType::IEEE_F32LE);
	dataset.close();
}

void Grid2DSphere::axis_vector(vector<Axis>& axes) const{
	axes = vector<Axis>({rAxis,thetaAxis});
}
double Grid2DSphere::zone_lorentz_factor(const int z_ind) const{
	abort(); // NOT IMPLEMENTED
}
void Grid2DSphere::get_connection_coefficients(EinsteinHelper* eh) const{ // default Minkowski
	PRINT_ASSERT(DO_GR,==,0);
	eh->christoffel.data = 0;
}
void Grid2DSphere::interpolate_shift(const double xup[4], double betaup[3], const unsigned dir_ind[NDIMS]) const{ // default Minkowski
	PRINT_ASSERT(DO_GR,==,0);
	betaup[0] = 0;
	betaup[1] = 0;
	betaup[2] = 0;
}
void Grid2DSphere::interpolate_3metric(const double xup[4], ThreeMetric* gammalow, const unsigned dir_ind[NDIMS]) const{ // default Minkowski
	PRINT_ASSERT(DO_GR,==,0);
	gammalow->data[ixx] = 1.0;
	gammalow->data[iyy] = 1.0;
	gammalow->data[izz] = 1.0;
	gammalow->data[ixy] = 0.0;
	gammalow->data[ixz] = 0.0;
	gammalow->data[iyz] = 0.0;
}
