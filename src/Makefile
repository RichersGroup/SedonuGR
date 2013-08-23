TARGET = ../gomc

OBJECTS = main.o Lua.o grid_general.o grid_1D_sphere.o grid_3D_cart.o emission.o scatter.o opacity.o locate_array.o cdf_array.o spectrum_array.o nulib_interface.o species_general.o photons.o transport.o neutrinos.o
TEST_OBJECTS = Lua.o grid_general.o grid_1D_sphere.o grid_3D_cart.o transport.o emission.o scatter.o opacity.o locate_array.o cdf_array.o spectrum_array.o nulib_interface.o species_general.o photons.o neutrinos.o

.PHONY: clean realclean

$(TARGET): $(OBJECTS) test
	$(CXX) $(CXXFLAGS) $(GSL_LDFLAGS) -o $(TARGET) $(OBJECTS) $(GSL_LIBS) $(NULIB_LIBS) $(HDF5LIBS)

.cc.o: 
	$(CXX) $(CXXFLAGS) $(GSL_INCL) -c -o -g $@ $<

.c.o: 
	cc -I. -c -o $@ $<

clean:
	rm -f $(OBJECTS)
	rm -f test.exe

realclean: clean
	rm -f $(TARGET)

test: $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(GSL_LDFLAGS) -o test.exe test.cpp $(TEST_OBJECTS) $(GSL_LIBS) $(NULIB_LIBS) $(HDF5LIBS)