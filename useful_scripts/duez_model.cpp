#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <cmath>

// Parameters
const int nx = 80;
const int ny = 80;
const int nz = 36;
const double xmin = -100.;
const double ymin = -100.;
const double zmin = 0;
const double xmax = 100.;
const double ymax = 100.;
const double zmax = 36.;
const int reflect_x = 0;
const int reflect_y = 0;
const int reflect_z = 1;

// Units
const double M_sun = 1.989e33;        // solar mass (g)
const double G     = 6.67259e-8;      // gravitational constant (cm^3/g/s^2)
const double c     = 2.99792458e10;   // speed of light (cm/s)
const double k_MeV = 8.6173324e-11;   // boltzmann constant (Mev/K)

using namespace std;

int main(int argc, char* argv[]){
  const int nvars = 6;
  if(argc != nvars+1){
    cout << "Usage: a.out rho T ye vx vy vz" << endl;
    exit(1);
  }

  // open the input/output files
  vector<ifstream> infiles(nvars);
  for(int i=0; i<infiles.size(); i++) infiles[i].open(argv[i+1]);
  ofstream out_mod;
  out_mod.open("model.mod");
  
  // write the metadata
  out_mod << "3D_cart" << endl;
  out_mod << "GRB" << endl;
  out_mod << nx << " " << ny << " " << nz << endl;
  out_mod << reflect_x << " " << reflect_y << " " << reflect_z << endl;
  out_mod << xmin << " " << xmax << endl;
  out_mod << ymin << " " << ymax << endl;
  out_mod << zmin << " " << zmax << endl;

  // first line is time - can ignore
  double trash;
  for(int i=0; i<infiles.size(); i++) infiles[i] >> trash;

  // read in the data
  vector<double> vals = vector<double>(nvars,0);
  for(int i=0; i<nx*ny*nz; i++){
    // read in the values
    for(int j=0; j<nvars; j++) infiles[j] >> vals[j];

    // convert units to CGS
    vals[0] *= pow(c,6)/ pow(G,3) / pow(M_sun,2);
    vals[1] /= k_MeV;
    vals[3] *= c;
    vals[4] *= c;
    vals[5] *= c;

    // write the values to the model file
    for(int j=0; j<nvars; j++) out_mod << vals[j] << "\t";
    out_mod << endl;
  }

  return 0;
}
