include make.inc
.PHONY: sedonu clean nulib nulibclean realclean

sedonu: #nulib
	mkdir -p exe
	$(MAKE) -C src
	ln -sf exe/sedonu

all: #nulib
	mkdir -p exe
	$(MAKE) all -C src
	ln -sf exe/sedonu

clean: 
	$(MAKE) -C src clean
	rm -rf exe

hdf5:
        wget https://support.hdfgroup.org/ftp/HDF5/current/src/hdf5-1.10.1.tar.gz
        tar -xvf hdf5-1.10.1.tar.gz
        rm -f hdf5-1.10.1.tar.gz
        mv hdf5-1.10.1 external/hdf5
        cd external/hdf5; ./configure --enable-fortran --enable-cxx; make; make install


realclean: clean nulibclean
	rm -rf exe

#########
# NuLib #
#########
nulib:
	$(MAKE) -C external NuLib

nulibclean:
	$(MAKE) -C external nulibclean
