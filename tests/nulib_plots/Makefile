.PHONY: all clean

rho=1e10  #g/cm^3
T=5      #MeV
Ye=0.3
E=5       #MeV
#species=1


all:
	mkdir -p species_0 species_1 species_2
	for species in 0 1 2 ; do \
	./nulib_test.exe ../../external/tables/NuLib.h5 $(rho) $(T) $(Ye) $(E) $$species ; \
	python eas_rho.py ; \
	python eas_T.py ; \
	python eas_ye.py ; \
	gnuplot eas_E_rho.gnuplot ; \
	gnuplot eas_E_T.gnuplot ; \
	gnuplot eas_E_ye.gnuplot ; \
	mv *.dat species_$$species ; \
	mv *.pdf species_$$species ; \
	done

clean:
	rm -rf species_0 species_1 species_2 plot_defaults.pyc *.dat