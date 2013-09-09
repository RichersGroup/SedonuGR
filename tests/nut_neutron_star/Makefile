.PHONY: all clean

all:
	python neutron_star.py > neutron_star.mod
	/usr/bin/mpirun -np 2 ./gomc
	gnuplot nut_neutron_star.gnuplot

clean:
	rm -f *.spec ray_* *.pdf neutron_star.mod