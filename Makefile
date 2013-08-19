include make.inc
.PHONY: clean

all: 
	cp make.inc external/NuLib
	$(MAKE) -C external/NuLib
	$(MAKE) -C src

clean: 
	$(MAKE) -C external/NuLib clean
	$(MAKE) -C src clean

realclean: clean
	$(MAKE) -C src realclean