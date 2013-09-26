.PHONY: all clean

all:
	python neutron_star.py > neutron_star.mod
	./gomc
	gnuplot equilibrium_ye.gnuplot

clean:
	rm -f *.spec ray_* *.pdf neutron_star.mod profile.dat