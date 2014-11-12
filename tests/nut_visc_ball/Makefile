.PHONY: all clean

all:
	python neutron_star.py > neutron_star.mod
	mpirun -np 1 ./gomc
	gnuplot compare_spectra.gnuplot

clean:
	rm -f spectrum_species*.dat ray_* fluid_* *.pdf neutron_star.mod
