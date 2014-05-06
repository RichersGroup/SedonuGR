.PHONY: all clean

all:
	rm -f output.dat
	cat neutron_star.py | sed 's/RHOHERE/1e+12/' | sed 's/TEMPHERE/1.0/' | sed 's/YEHERE/0.035/' > tmp.py
	python tmp.py > neutron_star.mod
	./gomc
	tail -n 1 fluid_00100 >> output.dat
	cat neutron_star.py | sed 's/RHOHERE/1e+12/' | sed 's/TEMPHERE/1.0/' | sed 's/YEHERE/0.1/' > tmp.py
	python tmp.py > neutron_star.mod
	./gomc
	tail -n 1 fluid_00100 >> output.dat
	cat neutron_star.py | sed 's/RHOHERE/1e+12/' | sed 's/TEMPHERE/1.0/' | sed 's/YEHERE/0.2/' > tmp.py
	python tmp.py > neutron_star.mod
	./gomc
	tail -n 1 fluid_00100 >> output.dat
	cat neutron_star.py | sed 's/RHOHERE/1e+12/' | sed 's/TEMPHERE/1.0/' | sed 's/YEHERE/0.3/' > tmp.py
	python tmp.py > neutron_star.mod
	./gomc
	tail -n 1 fluid_00100 >> output.dat
	cat neutron_star.py | sed 's/RHOHERE/1e+12/' | sed 's/TEMPHERE/1.0/' | sed 's/YEHERE/0.4/' > tmp.py
	python tmp.py > neutron_star.mod
	./gomc
	tail -n 1 fluid_00100 >> output.dat
	cat neutron_star.py | sed 's/RHOHERE/1e+12/' | sed 's/TEMPHERE/1.0/' | sed 's/YEHERE/0.5/' > tmp.py
	python tmp.py > neutron_star.mod
	./gomc
	tail -n 1 fluid_00100 >> output.dat
	cat neutron_star.py | sed 's/RHOHERE/1e+12/' | sed 's/TEMPHERE/1.0/' | sed 's/YEHERE/0.55/' > tmp.py
	python tmp.py > neutron_star.mod
	./gomc
	tail -n 1 fluid_00100 >> output.dat
	rm tmp.py

clean:
	rm -f *.spec ray_* fluid_* *.pdf neutron_star.mod tmp.py