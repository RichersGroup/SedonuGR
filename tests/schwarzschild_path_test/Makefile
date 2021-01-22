.PHONY: all clean

all:
	python3 empty_sphere.py 4 > empty_sphere.mod
	../../exe/schwarzschild_path_test radial1.lua > radial_smallstep4.dat
	../../exe/schwarzschild_path_test around1.lua > around_smallstep4.dat
	python3 check_results.py
	mv empty_sphere.mod empty_sphere_smallstep4.mod
	python3 empty_sphere.py 4 > empty_sphere.mod
	../../exe/schwarzschild_path_test radial2.lua > radial_smallstep2.dat
	../../exe/schwarzschild_path_test around2.lua > around_smallstep2.dat
	mv empty_sphere.mod empty_sphere_smallstep2.mod
	python3 empty_sphere.py 1 > empty_sphere.mod
	../../exe/schwarzschild_path_test radial1.lua > radial1.dat
	../../exe/schwarzschild_path_test around1.lua > around1.dat
	mv empty_sphere.mod empty_sphere1.mod
	python3 empty_sphere.py 2 > empty_sphere.mod
	../../exe/schwarzschild_path_test radial2.lua > radial2.dat
	../../exe/schwarzschild_path_test around2.lua > around2.dat
	mv empty_sphere.mod empty_sphere2.mod
	python3 empty_sphere.py 4 > empty_sphere.mod
	../../exe/schwarzschild_path_test radial4.lua > radial4.dat
	../../exe/schwarzschild_path_test around4.lua > around4.dat
	mv empty_sphere.mod empty_sphere4.mod
#	gnuplot makeplots.gpl
#	python makeplots.py

clean:
	rm -f spectrum_* ray_* fluid_* *.pdf *.dat empty_sphere.mod *~
