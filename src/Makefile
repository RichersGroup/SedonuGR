
CXX=mpicxx -O3 

GSL_BASE=$(GSL_DIR)
GSL_INCL=-I$(GSL_BASE)/include
GSL_LDFLAGS=-L$(GSL_BASE)/lib
GSL_LIBS=-lgsl -lgslcblas -llua

TARGET = gomc

OBJECTS = species_general.o photons.o main.o Lua.o grid_general.o grid_1D_sphere.o grid_3D_cart.o transport.o emission.o scatter.o opacity.o locate_array.o cdf_array.o spectrum_array.o

$(TARGET): $(OBJECTS)
	$(CXX) $(GSL_LDFLAGS) -o $(TARGET) $(OBJECTS) $(GSL_LIBS)

.cc.o: 
	$(CXX) $(GSL_INCL) -c -o -g $@ $<

.c.o: 
	cc -I. -c -o $@ $<

clean:
	rm -f $(OBJECTS) gomc

realclean: clean
	rm -f $(TARGET)
