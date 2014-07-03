#Only enable the commented compatability flags if compiling/running fails. Some systems (e.g. Hopper) won't compile with some of them enabled.

export
#general options
TARGET=$(CURDIR)/gomc
BASEDIR=$(CURDIR)

#compiler configuration
F90=mpif90
CXX=mpicxx
CC=mpicc

F90FLAGS= -O3 -fopenmp
CXXFLAGS= -O3 -fopenmp
CCFLAGS= -O3

#GSL
GSL_BASE=$(GSL_DIR)
GSL_INCL=-I$(GSL_BASE)/include
GSL_LIBS=-L$(GSL_BASE)/lib -lgsl -lgslcblas

#LUA
#LUA_BASE=$(CURDIR)/external/lua-5.2.2/install
LUA_INCL=-I/usr/include/lua5.2
LUA_LIBS=-llua5.2

#HDF5
HDF5_DIR=/usr/local/hdf5-1.8.13_gnu-4.8.2
HDF5INCS=-I$(HDF5_DIR)/include
HDF5LIBS=-L$(HDF5_DIR)/lib -lhdf5 -lhdf5_fortran -Wl,-rpath,$(HDF5_DIR)/lib

#NuLib
NULIB_DIR=$(CURDIR)/external/NuLib
NULIB_LIBS=$(NULIB_DIR)/src/nulib.a -lgfortran

#OTHER
OTHER_LIBS=