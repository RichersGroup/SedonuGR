.PHONY: all nulibclean NuLib

all: hdf5 NuLib

########
# HDF5 #
########
hdf5:
	wget http://www.hdfgroup.org/ftp/HDF5/current/src/hdf5-1.8.11.tar.gz
	tar -xvf hdf5-1.8.11.tar.gz
	cd hdf5-1.8.11; FC=$(F90) CC=$(CC) ./configure --enable-fortran
	cd hdf5-1.8.11; make
	cd hdf5-1.8.11; make test
	cd hdf5-1.8.11; make install
	mv hdf5-1.8.11/hdf5 .
	rm -rf hdf5-1.8.11*


#########
# NuLib #
#########
NuLib:
	echo > NuLib/make.inc # trick nulib into thinking make.inc exists there
	$(MAKE) -j1 -C NuLib

nulibclean:
	$(MAKE) -C NuLib clean
	rm -f NuLib/make.inc