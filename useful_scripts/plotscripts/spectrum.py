import matplotlib.pyplot as plt
import sys
import numpy as np
from pylab import *

##########
# INPUTS #
##########
dataFile = "/mnt/scratch/Christian_M1_Simulations/1D/spectrum_species2_00001.dat"

##############
# NOT INPUTS #
##############
h_MeV = 6.582119514e-22

# read the data
nu   = np.genfromtxt(dataFile,usecols=(0))
flux = np.genfromtxt(dataFile,usecols=(2))/8.

# set up the plot
plt.plot(h_MeV*nu,flux)
plt.axes().grid(True)
plt.axes().set_xlabel("Neutrino Energy (MeV)")
plt.axes().set_ylabel("Neutrino Flux (erg/s/Hz/sr)")
plt.show()
    
