.PHONY: all clean

all:
#	python neutron_star.py > neutron_star.mod
#	./gomc
	tail -n +2 neutron_star.mod > tmp.dat
	paste ray_00001 tmp.dat > profile.dat
	rm -f tmp.dat
	gnuplot equilibrium_ye.gnuplot

clean:
	rm -f *.spec ray_* *.pdf neutron_star.mod profile.dat