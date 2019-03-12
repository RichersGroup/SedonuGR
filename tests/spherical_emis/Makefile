.PHONY: all clean

all:
	python3 empty_sphere.py > empty_sphere.mod
	../../sedonu
	python3 plot_spectra.py

clean:
	rm -f spectrum_* ray_* fluid_* *.pdf empty_sphere.mod
