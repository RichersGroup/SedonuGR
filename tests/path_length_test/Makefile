.PHONY: all clean

#1 - absorption fraction
#2 - number of tau bins
#3 - max of tau bins
#4 - number of samples
#5 - number of scatters for each particle

all:
	../../exe/path_length_test 0.0 100 10 10000000 1 > output.dat
	gnuplot compare.gnuplot

clean:
	rm -f *.pdf output.dat
