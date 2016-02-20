include make.inc
.PHONY: all clean realclean tests testsclean hdf5 eos

all: 
	$(MAKE) -C src

clean: 
	$(MAKE) -C src clean

realclean: clean testsclean nulibclean eosclean
	$(MAKE) -C src realclean


#########
# TESTS #
#########

tests: all eos
	$(MAKE) tests -C src

testsclean:
	$(MAKE) testsclean -C src


#########
# NuLib #
#########
nulib:
	$(MAKE) -C external NuLib

nulibclean:
	$(MAKE) -C external nulibclean
