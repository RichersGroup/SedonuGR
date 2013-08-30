TARGET = ../gomc

HEADERS = cdf_array.h Error.h grid_1D_sphere.h grid_3D_cart.h grid_general.h
OBJECTS = main.o Lua.o grid_general.o grid_1D_sphere.o grid_3D_cart.o emission.o scatter.o opacity.o locate_array.o cdf_array.o spectrum_array.o nulib_interface.o species_general.o photons.o transport.o neutrinos.o
TEST_OBJECTS = Lua.o grid_general.o grid_1D_sphere.o grid_3D_cart.o transport.o emission.o scatter.o opacity.o locate_array.o cdf_array.o spectrum_array.o nulib_interface.o species_general.o photons.o neutrinos.o

.PHONY: clean realclean

# pull in dependency info for *existing* .o files so objects get rebuilt if header changes
-include $(OBJECTS:.o=.d)

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(GSL_LDFLAGS) -o $(TARGET) $(OBJECTS) $(GSL_LIBS) $(NULIB_LIBS) $(HDF5LIBS)

$(OBJECTS): %.o: %.cpp
	$(CXX) $(CXXFLAGS) $(GSL_INCL) -c -o $@ $<
	$(CXX) $(CXFLAGS) -MM $*.cpp > $*.d

clean:
	rm -f $(OBJECTS) test.exe *.d

realclean: clean
	rm -f $(TARGET)

test: $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(GSL_LDFLAGS) -o test.exe test.cpp $(TEST_OBJECTS) $(GSL_LIBS) $(NULIB_LIBS) $(HDF5LIBS)