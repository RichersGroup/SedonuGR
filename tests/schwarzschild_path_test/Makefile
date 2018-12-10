.PHONY: all clean

all:
	python3 empty_sphere.py > empty_sphere.mod
	../../exe/schwarzschild_path_test radial.lua > radial.dat
	../../exe/schwarzschild_path_test around.lua > around.dat
	python3 check_results.py
	gnuplot makeplots.gpl

clean:
	rm -f spectrum_* ray_* fluid_* *.pdf *.dat empty_sphere.mod *~
