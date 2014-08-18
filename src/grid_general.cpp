#pragma warning disable 161
#include <mpi.h>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
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
  ofstream outf;
  transport::open_file("fluid",iw,outf);
  outf << setprecision(4);
  outf << scientific;
  outf << "# r[0]\tr[1]\tr[2]\te_rad\trho\tT_gas\tYe\tt_therm\tt_lep" << endl;

  for (int i=0;i<z.size();i++)
  {
    double r[3];
    coordinates(i,r); 
    outf << r[0] << "\t";
    outf << r[1] << "\t";
    outf << r[2] << "\t";

    outf << z[i].e_rad << "\t";
    outf << z[i].rho << "\t";
    outf << z[i].T_gas*pc::k_MeV << "\t";
    outf << z[i].Ye << "\t";
    outf << 1.0 / fabs(1.0/z[i].t_eabs - 1.0/z[i].t_eemit) << "\t";
    outf << 1.0 / fabs(1.0/z[i].t_labs - 1.0/z[i].t_lemit) << "\t";
    outf << endl;
  }
  outf.close();
}
