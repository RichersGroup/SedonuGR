SOURCES=nuc_eos_press.cc testeos.cc readtable.cc nuc_eos_press_cs2.cc nuc_eos_dpdrhoe_dpderho.cc nuc_eos_short.cc nuc_eos_full.cc
CSOURCES= 
INCLUDES=nuc_eos.hh helpers.hh 
HDF5_DIR=/usr/local/hdf5-1.8.17_gnu-5.4.0/hdf5
HDF5INCS=-I$(HDF5_DIR)/include
HDF5LIBS=-L$(HDF5_DIR)/lib -lhdf5
GSL_BASE=$(GSL_DIR)
GSL_INCL=-I$(GSL_BASE)/include
GSL_LIBS=-L$(GSL_BASE)/lib -lgsl -lgslcblas
OBJECTS=$(SOURCES:.cc=.o )
COBJECTS=$(CSOURCES:.c=.o )
#EXTRAINCS=-Ivecmathlib

CC=gcc-4.8 
CPPC=g++-4.8

#CC=icc
#CPPC=icpc

#CC=gcc-mp-4.7 
#CPPC=g++-mp-4.7 
CPPFLAGS=-g -O0 -ansi -Wall -fopenmp -std=c++11 $(EXTRAINCS) -DHAVEGR
#CPPFLAGS=-g -ansi -Wall -ftree-vectorizer-verbose=6  -Winline -mtune=corei7 -fopenmp -msse4.2 -std=c++11 -Drestrict=__restrict__ -Ofast $(EXTRAINCS) -DHAVEGR
#CPPFLAGS=-g -ansi -Wall -ftree-vectorizer-verbose=6  -Winline -mtune=corei7 -fopenmp -msse2 -std=c++11 -Drestrict=__restrict__ -march=native -Ofast $(EXTRAINCS) -DHAVEGR
#-DDEBUG
CFLAGS=-g -O0 -std=gnu99 -ansi -Wall -fopenmp
EXTRALIBS=-lm -L/usr/local/lib -lgsl -lgslcblas -Wl,-rpath,$(HDF5DIR)/lib

all: nuc_eos.a Munue_of_RhoTempYe Ye_of_RhoTempMunue All_of_RhoTempYe many_Ye

many_Ye: $(OBJECTS) $(COBJECTS) $(INCLUDES)
	$(CPPC) $(CPPFLAGS) $(HDF5INCS) -c many_Ye.cc -o many_Ye.o
	$(CPPC) $(CPPFLAGS) -o many_Ye $(OBJECTS) $(COBJECTS) $(HDF5LIBS) $(EXTRALIBS) many_Ye.o

Munue_of_RhoTempYe: $(OBJECTS) $(COBJECTS) $(INCLUDES)
	$(CPPC) $(CPPFLAGS) $(HDF5INCS) -c Munue_of_RhoTempYe.cc -o Munue_of_RhoTempYe.o
	$(CPPC) $(CPPFLAGS) -o Munue_of_RhoTempYe $(OBJECTS) $(COBJECTS) $(HDF5LIBS) $(EXTRALIBS) Munue_of_RhoTempYe.o

Ye_of_RhoTempMunue: $(OBJECTS) $(COBJECTS) $(INCLUDES)
	$(CPPC) $(CPPFLAGS) $(HDF5INCS) -c Ye_of_RhoTempMunue.cc -o Ye_of_RhoTempMunue.o
	$(CPPC) $(CPPFLAGS) -o Ye_of_RhoTempMunue $(OBJECTS) $(COBJECTS) $(HDF5LIBS) $(EXTRALIBS) Ye_of_RhoTempMunue.o

All_of_RhoTempYe: $(OBJECTS) $(COBJECTS) $(INCLUDES)
	$(CPPC) $(CPPFLAGS) $(HDF5INCS) -c All_of_RhoTempYe.cc -o All_of_RhoTempYe.o
	$(CPPC) $(CPPFLAGS) -o All_of_RhoTempYe $(OBJECTS) $(COBJECTS) $(HDF5LIBS) $(EXTRALIBS) All_of_RhoTempYe.o

testi: testi.cc
	$(CPPC) $(CPPFLAGS) -o testi testi.cc -Ivecmathlib

$(OBJECTS): %.o: %.cc $(INCLUDES)
	$(CPPC) $(CPPFLAGS) $(HDF5INCS) -c $< -o $@


$(COBJECTS): %.o: %.c $(INCLUDES)
	$(CC) $(CFLAGS) $(HDF5INCS) -c $< -o $@

nuc_eos.a: $(OBJECTS)
	ar r nuc_eos.a nuc_eos*.o readtable.o

clean:
	rm -f *.o
	rm -rf *.cc.*
	rm Munue_of_RhoTempYe Ye_of_RhoTempMunue All_of_RhoTempYe many_Ye
