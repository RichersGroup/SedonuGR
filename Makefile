include make.inc
.PHONY: sedonu clean nulib nulibclean realclean

sedonu: #nulib
	mkdir -p exe
	$(MAKE) -C src
	ln -sf exe/sedonu

clean: 
	$(MAKE) -C src clean

realclean: clean nulibclean
	rm -rf exe

#########
# NuLib #
#########
nulib:
	$(MAKE) -C external NuLib

nulibclean:
	$(MAKE) -C external nulibclean
