#include <mpi.h>
#include <math.h>
#include <iostream>
#include <iomanip>
#include "grid_general.h"
#include "physical_constants.h"
#include <fstream>
#include <gsl/gsl_rng.h>
#include "Lua.h"

namespace pc = physical_constants;
using namespace std;


//------------------------------------------------------------
// initialize the grid
//------------------------------------------------------------
void grid_general::init(Lua* lua)
{
  cout << "in grid_general::init" << endl;
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
// Combine the radiation tallies in all zones
// from all processors using MPI allreduce
//------------------------------------------------------------
void grid_general::reduce_radiation()
{
  // largest block to reduce
  int bsize    = 10000;
  int n_blocks = floor(z.size()/bsize);

  // reduce in blocks
  for (int i=0;i<n_blocks;i++) 
    reduce_radiation_block(bsize,i*bsize);
  
  // get remainder
  int remainder = z.size() - n_blocks*bsize;
  if (remainder > 0) reduce_radiation_block(remainder,n_blocks*bsize);
}


//------------------------------------------------------------
// Combine the radiation tallies in blocks
// from all processors  using MPI allreduce
//------------------------------------------------------------
void grid_general::reduce_radiation_block(int bsize, int start)
{
  int j,size;
  MPI_Comm_size( MPI_COMM_WORLD, &size );

  double *src = new double[bsize];
  double *dst = new double[bsize];
  
  // reduce (average) e_rad
  for (j=0;j<bsize;j++) 
  {
    src[j] = z[start + j].e_rad/size;
    dst[j] = 0;
  }
  MPI_Allreduce(src,dst,bsize,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  for (j=0;j<bsize;j++) z[start + j].e_rad = dst[j];


  // reduce (average) e_abs
  for (j=0;j<bsize;j++) 
  {
    src[j] = z[start + j].e_abs/size;
    dst[j] = 0;
  }
  MPI_Allreduce(src,dst,bsize,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  for (j=0;j<bsize;j++) z[start + j].e_abs = dst[j];

  // TODO - need to put in other quantities...
  // TODO - perhaps put the radiation arrays in transport?
  delete src;
  delete dst;

}


void grid_general::reduce_T_gas()
{
  // reduce gas temperature 
  double *src_ptr = new double[z.size()];
  double *dst_ptr = new double[z.size()];

  for (int i=0;i<z.size();i++) {src_ptr[i] = z[i].T_gas; dst_ptr[i] = 0.0;}
  MPI_Allreduce(src_ptr,dst_ptr,z.size(),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  for (int i=0;i<z.size();i++) z[i].T_gas = dst_ptr[i];

  delete src_ptr;
  delete dst_ptr;

}


void grid_general::write_zones(int iw)
{
  char zonefile[1000];
  char base[1000];

  if (iw < 10) sprintf(base,"_0000%d",iw);
  else if (iw < 100) sprintf(base,"_000%d",iw);
  else if (iw < 1000) sprintf(base,"_00%d",iw);
  else if (iw < 10000) sprintf(base,"_0%d",iw);
  else sprintf(base,"_%d",iw);
  sprintf(zonefile,"ray%s",base);

  ofstream outf;
  outf.open(zonefile);
  //  outf << setw(12);
  outf << setprecision(4);
  outf << scientific;

  for (int i=0;i<z.size();i++)
  {
    double r[3];
    coordinates(i,r); 
    outf << r[0] << " ";

    double T_rad = pow(z[i].e_rad/pc::a,0.25);
    outf << T_rad << " ";
    outf << z[i].T_gas << " ";
    
    outf << endl;
  }

}
