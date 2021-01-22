all:
	python3 uniform_sphere.py > uniform_sphere.mod
	../../sedonu param.lua
	python3 uniform_sphere_test.py

clean:
	rm -f fluid_*.h5 oven.mod *.pdf
