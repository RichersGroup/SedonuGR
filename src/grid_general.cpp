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

  // set the reduction block size
  block_size  = lua->scalar<int>("block_size");

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
  int n_blocks = floor(z.size()/block_size);

  // reduce in blocks
  for (int i=0;i<n_blocks;i++) 
    reduce_radiation_block(block_size,i*block_size);
  
  // get remainder
  int remainder = z.size() - n_blocks*block_size;
  if (remainder > 0) reduce_radiation_block(remainder,n_blocks*block_size);
}


//------------------------------------------------------------
// Combine the radiation tallies in blocks
// from all processors  using MPI allreduce
//------------------------------------------------------------
void grid_general::reduce_radiation_block(int block_size, int start)
{
  int nprocs;
  MPI_Comm_size( MPI_COMM_WORLD, &nprocs );

  vector<real> send(block_size,0);
  vector<real> receive(block_size,0);
  MPI_Datatype type = ( sizeof(real)==4 ? MPI_FLOAT : MPI_DOUBLE);

  // average e_rad
  #pragma omp parallel for
  for(int i=0; i<block_size; i++) send[i] = z[start+i].e_rad/nprocs;
  MPI_Allreduce(&send.front(), &receive.front(), block_size, type, MPI_SUM, MPI_COMM_WORLD);
  #pragma omp parallel for
  for(int i=0; i<block_size; i++) z[start+i].e_rad = receive[i];

  // average e_abs
  #pragma omp parallel for
  for(int i=0; i<block_size; i++) send[i] = z[start+i].e_abs/nprocs;
  MPI_Allreduce(&send.front(), &receive.front(), block_size, type, MPI_SUM, MPI_COMM_WORLD);
  #pragma omp parallel for
  for(int i=0; i<block_size; i++) z[start+i].e_rad = receive[i];

  // average l_abs
  #pragma omp parallel for
  for(int i=0; i<block_size; i++) send[i] = z[start+i].l_abs/nprocs;
  MPI_Allreduce(&send.front(), &receive.front(), block_size, type, MPI_SUM, MPI_COMM_WORLD);
  #pragma omp parallel for
  for(int i=0; i<block_size; i++) z[start+i].e_rad = receive[i];

  // TODO - need to put in other quantities...
  // TODO - can optimize by having only one for loop before and one loop after, but uses more memory
}


void grid_general::reduce_Ye()
{
  vector<real> send(z.size(),0);
  vector<real> receive(z.size(),0);
  MPI_Datatype type = ( sizeof(real)==4 ? MPI_FLOAT : MPI_DOUBLE);

  // reduce Ye
  #pragma omp parallel for
  for (int i=0;i<z.size();i++) send[i] = z[i].Ye;
  MPI_Allreduce(&send.front(), &receive.front(), z.size(), type, MPI_SUM, MPI_COMM_WORLD);
  #pragma omp parallel for
  for (int i=0;i<z.size();i++) z[i].Ye = receive[i];
}

void grid_general::reduce_T()
{
  vector<real> send(z.size(),0);
  vector<real> receive(z.size(),0);
  MPI_Datatype type = ( sizeof(real)==4 ? MPI_FLOAT : MPI_DOUBLE);

  // reduce gas temperature
  #pragma omp parallel for
  for (int i=0;i<z.size();i++) send[i] = z[i].T_gas;
  MPI_Allreduce(&send.front(), &receive.front(), z.size(), type, MPI_SUM, MPI_COMM_WORLD);
  #pragma omp parallel for
  for (int i=0;i<z.size();i++) z[i].T_gas = receive[i];
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
    outf << z[i].Ye << " ";
    outf << endl;
  }

}
