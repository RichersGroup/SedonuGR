.PHONY: all clean

all: gomc
	python make_mod.py > lucy.mod
	./gomc
	gnuplot compare_spectrum.gnuplot

clean:
	rm -f *.spec ray_* lucy.mod compare_spectrum.pdf

realclean: clean
	rm -f gomc