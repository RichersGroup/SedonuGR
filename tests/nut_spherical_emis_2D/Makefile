.PHONY: all clean

all:
	python empty_sphere.py > neutron_star.mod
	../../sedonu
	gnuplot compare_spectra.gnuplot

clean:
	rm -f spectrum_* ray_* fluid_* *.pdf empty_sphere.mod
