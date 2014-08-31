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
	vector<double> r;
	zone_coordinates(0,r);
	int dimensionality = r.size();

	ofstream outf;
	transport::open_file("fluid",iw,outf);
	zone::write_header(dimensionality,outf);
	vector<int> dir_ind;

	for (unsigned z_ind=0; z_ind<z.size(); z_ind++)
	{
	        zone_directional_indices(z_ind, dir_ind);
		if(dir_ind[dir_ind.size()-1]==0) outf << endl;
		zone_coordinates(z_ind,r);
		z[z_ind].write_line(r,outf);
	}
	outf.close();
}
