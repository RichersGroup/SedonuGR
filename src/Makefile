TARGET = ../gomc

OBJECTS = main.o Lua.o grid_general.o grid_1D_sphere.o grid_3D_cart.o emission.o scatter.o opacity.o locate_array.o cdf_array.o spectrum_array.o nulib_interface.o species_general.o photons.o transport.o neutrinos.o
TEST_OBJECTS = Lua.o grid_general.o grid_1D_sphere.o grid_3D_cart.o transport.o emission.o scatter.o opacity.o locate_array.o cdf_array.o spectrum_array.o nulib_interface.o species_general.o photons.o neutrinos.o

.PHONY: clean realclean

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(GSL_LDFLAGS) -o $(TARGET) $(OBJECTS) $(GSL_LIBS) $(NULIB_LIBS) $(HDF5LIBS)

$(OBJECTS): %.o: %.cpp
	$(CXX) $(CXXFLAGS) $(GSL_INCL) -c -o $@ $<
	$(CXX) $(CXFLAGS) -MM $*.cpp > $*.d

# pull in dependency info for *existing* .o files so objects get rebuilt if header changes
-include $(OBJECTS:.o=.d)

clean:
	rm -f $(OBJECTS) test.exe *.d

realclean: clean
	rm -f $(TARGET)

tests: $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(GSL_LDFLAGS) -o ../tests/nulib_plots/nulib_test.exe nulib_test.cpp $(TEST_OBJECTS) $(GSL_LIBS) $(NULIB_LIBS) $(HDF5LIBS)
	cp $(TARGET) ../tests/phot_spherical_emis
	cp $(TARGET) ../tests/nut_spherical_emis
	cp $(TARGET) ../tests/nut_neutron_star