TARGET = ../gomc

OBJECTS = main.o Lua.o grid_general.o grid_0D_isotropic.o grid_1D_sphere.o grid_2D_sphere.o grid_3D_cart.o emission.o scatter.o opacity.o locate_array.o cdf_array.o spectrum_array.o nulib_interface.o species_general.o photons.o transport.o neutrinos.o thread_RNG.o solve_equilibrium.o propagate.o relativity.o zone.o particle.o tools.o
TEST_OBJECTS = Lua.o grid_general.o grid_0D_isotropic.o grid_1D_sphere.o grid_2D_sphere.o grid_3D_cart.o transport.o emission.o scatter.o opacity.o locate_array.o cdf_array.o spectrum_array.o nulib_interface.o species_general.o photons.o neutrinos.o thread_RNG.o solve_equilibrium.o propagate.o relativity.o zone.o particle.o tools.o

.PHONY: clean realclean

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJECTS) $(GSL_LIBS) $(LUA_LDFLAGS) $(NULIB_LIBS) $(HDF5LIBS) $(LUA_LIBS) $(OTHER_LIBS)

$(OBJECTS): %.o: %.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<  $(GSL_INCL) $(LUA_INCL) $(OTHER_INCL) $(HDF5INCS)

# pull in dependency info for *existing* .o files so objects get rebuilt if header changes
-include $(OBJECTS:.o=.d)

clean:
	rm -f $(OBJECTS) test.exe *.d

realclean: clean
	rm -f $(TARGET)

tests: $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o ../tests/nulib_plots/nulib_test.exe nulib_test.cpp $(TEST_OBJECTS) $(GSL_LIBS) $(NULIB_LIBS) $(HDF5LIBS) $(LUA_LIBS) $(OTHER_LIBS) $(LUA_INCL)
	$(CXX) $(CXXFLAGS) -o ../tests/nulib_plots/nulib_eas_single nulib_eas_single.cpp $(TEST_OBJECTS) $(GSL_LIBS) $(NULIB_LIBS) $(HDF5LIBS) $(LUA_LIBS) $(OTHER_LIBS) $(LUA_INCL)
	$(CXX) $(CXXFLAGS) -o ../tests/nut_blackbody/nut_blackbody.exe nut_blackbody.cpp $(TEST_OBJECTS) $(GSL_LIBS) $(NULIB_LIBS) $(HDF5LIBS) $(LUA_LIBS) $(OTHER_LIBS) $(LUA_INCL)
	$(CXX) $(CXXFLAGS) -o ../tests/nut_therm_shell/nut_therm_shell.exe nut_therm_shell.cpp $(TEST_OBJECTS) $(GSL_LIBS) $(NULIB_LIBS) $(HDF5LIBS) $(LUA_LIBS) $(OTHER_LIBS) $(LUA_INCL)
	cp $(TARGET) ../tests/phot_spherical_emis
	cp $(TARGET) ../tests/phot_visc_ball
	cp $(TARGET) ../tests/nut_spherical_emis
	cp $(TARGET) ../tests/nut_spherical_emis_2D
	cp $(TARGET) ../tests/nut_visc_ball
	cp $(TARGET) ../tests/nut_neutron_star
	cp $(TARGET) ../tests/nut_neutron_star_2D
	cp $(TARGET) ../tests/equilibrium_ye

testsclean:
	$(MAKE) -C ../tests/phot_spherical_emis clean
	$(MAKE) -C ../tests/phot_visc_ball clean
	$(MAKE) -C ../tests/nut_spherical_emis_2D clean
	$(MAKE) -C ../tests/nut_spherical_emis clean
	$(MAKE) -C ../tests/nut_visc_ball clean
	$(MAKE) -C ../tests/nut_therm_shell clean
	$(MAKE) -C ../tests/nut_neutron_star clean
	$(MAKE) -C ../tests/nut_neutron_star_2D clean
	$(MAKE) -C ../tests/nut_blackbody clean
	$(MAKE) -C ../tests/equilibrium_ye clean
	$(MAKE) -C ../tests/nulib_plots clean
	rm -f ../tests/*/gomc
	rm -f ../tests/nulib_plots/nulib_test.exe
