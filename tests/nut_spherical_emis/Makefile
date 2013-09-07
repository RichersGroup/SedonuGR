.PHONY: all clean

all:
	python empty_sphere.py > empty_sphere.mod
	./gomc
	gnuplot compare_spectra.gnuplot

clean:
	rm -f *.spec ray_* *.pdf empty_sphere.mod