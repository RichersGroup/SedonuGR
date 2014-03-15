SOURCES=driver.cc nuc_eos_press.cc testeos.cc readtable.cc nuc_eos_press_cs2.cc nuc_eos_dpdrhoe_dpderho.cc nuc_eos_short.cc nuc_eos_full.cc
CSOURCES= 
INCLUDES=nuc_eos.hh helpers.hh 
HDF5DIR=/usr/local/hdf5-1.8.11_gnu-4.6.3
HDF5INCS=-I$(HDF5DIR)/include
HDF5LIBS=-L$(HDF5DIR)/lib -lhdf5
GSL_BASE=$(GSL_DIR)
GSL_INCL=-I$(GSL_BASE)/include
GSL_LIBS=-L$(GSL_BASE)/lib -lgsl -lgslcblas
OBJECTS=$(SOURCES:.cc=.o )
COBJECTS=$(CSOURCES:.c=.o )
#EXTRAINCS=-Ivecmathlib

#CC=gcc-4.8 
#CPPC=g++-4.8

CC=icc
CPPC=icpc

#CC=gcc-mp-4.7 
#CPPC=g++-mp-4.7 
CPPFLAGS=-g -O0 -ansi -Wall -fopenmp -std=c++11 $(EXTRAINCS) -DHAVEGR
#CPPFLAGS=-g -ansi -Wall -ftree-vectorizer-verbose=6  -Winline -mtune=corei7 -fopenmp -msse4.2 -std=c++11 -Drestrict=__restrict__ -Ofast $(EXTRAINCS) -DHAVEGR
#CPPFLAGS=-g -ansi -Wall -ftree-vectorizer-verbose=6  -Winline -mtune=corei7 -fopenmp -msse2 -std=c++11 -Drestrict=__restrict__ -march=native -Ofast $(EXTRAINCS) -DHAVEGR
#-DDEBUG
CFLAGS=-g -O0 -std=gnu99 -ansi -Wall -fopenmp
EXTRALIBS=-lm -L/usr/local/lib -lgsl -lgslcblas -Wl,-rpath,$(HDF5DIR)/lib

all: driver nuc_eos.a

testi: testi.cc
	$(CPPC) $(CPPFLAGS) -o testi testi.cc -Ivecmathlib

driver: $(OBJECTS) $(COBJECTS) $(INCLUDES)
	$(CPPC) $(CPPFLAGS) -o driver $(OBJECTS) $(COBJECTS) $(HDF5LIBS) $(EXTRALIBS)

$(OBJECTS): %.o: %.cc $(INCLUDES)
	$(CPPC) $(CPPFLAGS) $(HDF5INCS) -c $< -o $@


$(COBJECTS): %.o: %.c $(INCLUDES)
	$(CC) $(CFLAGS) $(HDF5INCS) -c $< -o $@

nuc_eos.a: $(OBJECTS)
	ar r nuc_eos.a nuc_eos*.o readtable.o

clean:
	rm -f *.o
	rm -rf *.cc.*
