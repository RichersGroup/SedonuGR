.PHONY: all clean

all:
	rm -f output.dat
	cat neutron_star.py | sed 's/RHOHERE/1e+12/' | sed 's/TEMPHERE/10/' | sed 's/YEHERE/0.3/' > tmp.py
	python tmp.py > neutron_star.mod
	./gomc
	rm tmp.py

clean:
	rm -f *.spec ray_* fluid_* *.pdf neutron_star.mod tmp.py