.PHONY: all clean

rho=1e12  #g/cm^3
T=10      #MeV
Ye=0.2
E=10      #MeV
mu=4.44614 #MeV


all:
	mkdir -p species_0 species_1 species_2
	for species in 0 1 2 ; do \
	./nulib_test.exe ../../external/tables/NuLib_LS220_rho150_temp90_ye60_ng24_ns3_version1.0_20140327.h5 $(rho) $(T) $(Ye) $(mu) $(E) $$species ; \
	python eas_rho.py ; \
	python eas_T.py ; \
	python eas_ye.py ; \
	gnuplot eas_E_rho.gnuplot ; \
	gnuplot eas_E_T.gnuplot ; \
	gnuplot eas_E_ye.gnuplot ; \
	gnuplot emis_vs_opac.gnuplot; \
	mv *.dat species_$$species ; \
	mv *.pdf species_$$species ; \
	done

clean:
	rm -rf species_0 species_1 species_2 plot_defaults.pyc *.dat