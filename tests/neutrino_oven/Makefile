all:
	echo "specify N or GR"

N:
	python3 oven_N.py > oven.mod
	../../sedonu param_abs.lua
	python3 oven_test.py
	../../sedonu param_mildscatter.lua
	python3 oven_test.py
	../../sedonu param_heavyscatter.lua
	python3 oven_test.py

GR:
	python3 oven_GR.py > oven.mod
	../../sedonu param_abs.lua
	python3 oven_test.py
	../../sedonu param_mildscatter.lua
	python3 oven_test.py
	../../sedonu param_heavyscatter.lua
	python3 oven_test.py

clean:
	rm -f fluid_*.h5 oven.mod *.pdf
