include make.inc
.PHONY: all clean realclean tests testsclean hdf5

all: 
	$(MAKE) -C src

clean: 
	$(MAKE) -C src clean

realclean: clean
	$(MAKE) -C src realclean


#########
# TESTS #
#########

tests: all
	$(MAKE) tests -C src

testsclean:
	$(MAKE) testsclean -C src


########
# HDF5 #
########
hdf5:
	$(MAKE) -C external hdf5


#########
# NuLib #
#########
nulib:
	$(MAKE) -C external NuLib

nulibclean:
	$(MAKE) -C external nulibclean
