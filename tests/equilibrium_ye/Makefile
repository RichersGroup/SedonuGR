.PHONY: all clean

all:
	python neutron_star.py > neutron_star.mod
	./gomc
#	../../external/EOSsuper/driver neutron_star.mod > eos.dat
#	cat eos.dat | tail -n +5 eos.dat > tmp.dat && mv tmp.dat eos.dat
	gnuplot equilibrium_ye.gnuplot

clean:
	rm -f spectrum_species*.dat ray_* fluid_* *.pdf neutron_star.mod profile.dat eos.dat
