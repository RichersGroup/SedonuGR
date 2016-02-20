.PHONY: all clean

rho=1e10  #g/cm^3
T=10      #MeV
Ye=0.1
E=10      #MeV
EOS="/data/tables/EOS/LS220.h5"
NULIB="/data/tables/NuLib/NuLib_LS220_fiducial.h5"

all:
	$(eval mu=$(shell ../../external/EOSsuper/Munue_of_RhoTempYe $(EOS) $(rho) $(T) $(Ye)))
	mkdir -p species_0 species_1 species_2
	for species in 0 1 2 ; do \
	../../exe/nulib_test $(NULIB) $(rho) $(T) $(Ye) $(mu) $(E) $$species ; \
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
