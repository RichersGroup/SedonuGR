#include <iostream>
#include <cstring>
#include <cstdlib>
#include "hdf5.h"
//#include "H5Cpp.h"
#include "hdf5_tools.h"

int main(){

  //check that the file exists
  const char* filename = "/data/MCNuTrans/Fernandez/torus_hdf5_plt_cnt_0280";
  std::cout << filename << std::endl;
  if (!file_is_readable(filename)) {
    fprintf(stderr, "Could not read table");
  }

  hid_t file;
  HDF5_ERROR(file = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT));

  hid_t dataset,dapl;
  HDF5_ERROR(dataset = H5Dopen(file, "efrc",dapl));
  DataSpace filespace = dataset.getSpace();
  HDF5_ERROR(H5Dclose(dataset));  

  HDF5_ERROR(H5Fclose(file));
  return 0;
}
