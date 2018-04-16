.PHONY: all clean

all:
	python neutron_star.py > neutron_star.mod
	../../sedonu
	../../useful_scripts/EOSsuper/many_Ye SFHo.h5 neutron_star.mod > eos.dat
	python plot.py

clean:
	rm -f spectrum_species*.dat ray_* fluid_* *.pdf neutron_star.mod profile.dat eos.dat
