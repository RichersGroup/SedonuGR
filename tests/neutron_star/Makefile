.PHONY: all clean

all:
	python neutron_star.py > neutron_star.mod
	../../sedonu
	python plot.py

clean:
	rm -f spectrum*.dat ray_* fluid_* *.pdf neutron_star.mod
