#pragma warning disable 161
#include <mpi.h>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
#include <cassert>
#include "grid_general.h"
#include "physical_constants.h"
#include "Lua.h"

namespace pc = physical_constants;
using namespace std;


//------------------------------------------------------------
// initialize the grid
//------------------------------------------------------------
void grid_general::init(Lua* lua)
{
  // read the model file or fill in custom model
  std::string model_file = lua->scalar<std::string>("model_file");
  if(model_file == "custom_model") custom_model(lua);
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

  for (int i=0;i<z.size();i++)
  {
    zone_coordinates(i,r); 
    z[i].write_line(r,outf);
  }
  outf.close();
}
