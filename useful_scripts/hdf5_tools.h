#ifndef _HDF5_TOOLS_H
#define _HDF5_TOOLS_H
#include <cstdio>

// Catch HDF5 errors
#define HDF5_ERROR(fn_call)				\
  do {							\
    int _error_code = fn_call;				\
    if (_error_code < 0) {				\
      fprintf(stderr,					\
	      "HDF5 call '%s' returned error code %d",	\
	      #fn_call, _error_code);			\
      abort(); }					\
  } while (0)
  
static int file_is_readable(const char* filename)	\
{							\
  FILE* fp = NULL;					\
  fp = fopen(filename, "r");				\
  if(fp != NULL)					\
    {							\
      fclose(fp);					\
      return 1;						\
    }							\
  return 0;						\
}							\

    
// Use these two defines to easily read in a lot of variables in the same way
// The first reads in one variable of a given type completely
#define READ_EOS_HDF5(NAME,VAR,TYPE,MEM)				\
  do {									\
    hid_t dataset;							\
    HDF5_ERROR(dataset = H5Dopen(file, NAME));				\
    HDF5_ERROR(H5Dread(dataset, TYPE, MEM, H5S_ALL, H5P_DEFAULT, VAR));	\
    HDF5_ERROR(H5Dclose(dataset));					\
  } while (0)
// The second reads a given variable into a hyperslab of the alltables_temp array
#define READ_EOSTABLE_HDF5(NAME,OFF)					\
  do {									\
    hsize_t offset[2]     = {OFF,0};					\
    H5Sselect_hyperslab(mem3, H5S_SELECT_SET, offset, NULL, var3, NULL); \
    READ_EOS_HDF5(NAME,alltables_temp,H5T_NATIVE_DOUBLE,mem3);		\
  } while (0)								\
    

#endif
