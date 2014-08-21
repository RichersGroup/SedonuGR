.PHONY: all clean

all:
	python empty_sphere.py > empty_sphere.mod
	./gomc
	gnuplot compare_spectra.gnuplot

clean:
	rm -f spectrum_* ray_* fluid_* *.pdf empty_sphere.mod
