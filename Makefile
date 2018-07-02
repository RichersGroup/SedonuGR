include make.inc
.PHONY: sedonu clean nulib nulibclean realclean lua

sedonu: #nulib
	mkdir -p exe
	$(MAKE) -C src
	ln -sf exe/sedonu

all: #nulib
	mkdir -p exe
	$(MAKE) all -C src
	ln -sf exe/sedonu

clean: 
	$(MAKE) -C src clean
	rm -rf exe

hdf5:
	if [ -d "external/hdf5/bin" ]; then \
		echo "HDF5 already exists"; \
	else \
		wget https://support.hdfgroup.org/ftp/HDF5/current/src/hdf5-1.10.1.tar.gz; \
		tar -xvf hdf5-1.10.1.tar.gz; \
		rm -f hdf5-1.10.1.tar.gz; \
		cd hdf5-1.10.1; FC=$(F90) ./configure --enable-fortran --enable-cxx; make -j; make install; cd ..; \
		mv hdf5-1.10.1/hdf5 external/; \
		rm -rf hdf5-1.10.1; \
	fi
lua:
	if [ -d "external/lua/bin" ]; then \
		echo "Lua already exists"; \
	else \
		wget http://www.lua.org/ftp/lua-5.3.4.tar.gz; \
		tar -xvf lua-5.3.4.tar.gz; \
		rm -f lua-5.3.4.tar.gz; \
		cd lua-5.3.4; make linux; make local; cd ..; \
		rm -rf external/lua; \
		mv lua-5.3.4/install external/lua; \
		rm -rf lua-5.3.4; \
	fi

realclean:
	$(MAKE) clean
	$(MAKE) nulibclean
	rm -rf exe external/hdf5

#########
# NuLib #
#########
nulib:
	$(MAKE) -C external NuLib

nulibclean:
	$(MAKE) -C external nulibclean
