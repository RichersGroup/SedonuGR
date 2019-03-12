all:
	python3 oven.py > oven.mod
	../../sedonu param.lua
	python3 oven_test.py
	../../sedonu param_abs.lua
	python3 oven_test.py
	../../sedonu param_mildscatter.lua
	python3 oven_test.py
	../../sedonu param_heavyscatter.lua
	python3 oven_test.py

clean:
	rm -f fluid_*.h5 oven.mod
