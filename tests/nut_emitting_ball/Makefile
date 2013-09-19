.PHONY: all clean

all:
	python neutron_star.py > neutron_star.mod
	./gomc
	gnuplot compare_spectra.gnuplot

clean:
	rm -f *.spec ray_* *.pdf neutron_star.mod