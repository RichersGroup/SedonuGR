.PHONY: clean

rho=1e14  #g/cm^3                                                                                                                      
T=1       #MeV                                                                                                                         
Ye=0.35
E=1       #MeV                                                                                                                         
species=0


all:
	./nulib_test.exe ../../external/tables/NuLib_LS220.h5 $(rho) $(T) $(Ye) $(E) $(species)	
	python eas_rho.py
	python eas_T.py
	python eas_ye.py
	gnuplot eas_E_rho.gnuplot
	gnuplot eas_E_T.gnuplot
	gnuplot eas_E_ye.gnuplot

clean:
	rm -f *.dat *.pdf