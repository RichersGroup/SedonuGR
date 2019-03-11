all:
	python3 oven.py > oven.mod
	../../sedonu param.lua > output.txt
	python3 oven_test.py
	../../sedonu param_abs.lua > output.txt
	python3 oven_test.py

clean:
	rm -f fluid_*.h5 oven.mod
