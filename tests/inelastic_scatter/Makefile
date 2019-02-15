.PHONY: all clean

all:
	../../exe/scattering_test elastic_isotropic_kernel.lua > elastic_isotropic_kernel.dat
	python3 check_results.py

clean:
	rm -f spectrum_* ray_* fluid_* *.pdf *.dat empty_sphere.mod *~
