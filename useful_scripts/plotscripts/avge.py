import matplotlib.pyplot as plt
import h5py
import sys
import numpy as np
from pylab import *

##########
# INPUTS #
##########
dataFile = "/mnt/scratch/Christian_M1_Simulations/3D/fluid_00001.h5"
species = 2 # 0:nue 1:nuebar 2:nux
group = -1 # -1:integrated
normalDirection = 2 # which direction is normal to the plotting plane? 0:x 1:y 2:z
normalIndex = 25 # index of zone in normalDirection direction

##############
# NOT INPUTS #
##############
h_MeV = 6.582119514e-22

# open the hdf5 file
f = h5py.File(dataFile,"r")

# read the grid structure and erad data
xgrid = np.array(f["grid_0(cm)"])/1e5
ygrid = np.array(f["grid_1(cm)"])/1e5
zgrid = np.array(f["grid_2(cm)"])/1e5
Egrid = np.array(f["distribution(erg|ccm,lab)/distribution_frequency_grid(Hz,lab)"])*h_MeV
invEmid = np.array([1.0/(Egrid[i]+Egrid[i+1]) for i in range(len(Egrid)-1) ])
M0 = np.array(f["distribution(erg|ccm,lab)/M0"])[:,:,:,species,:,0]
N0 = M0 * invEmid

plotdata = np.sum(M0,axis=3) / np.sum(N0,axis=3)
plotdata[np.where(plotdata!=plotdata)] = 1
plotdata = log10(plotdata)
print(np.min(plotdata),np.max(plotdata))

# select out only the plane we want
if normalDirection==0:
    plotdata_2D = plotdata[normalIndex,:,:]
    grid1 = ygrid
    grid2 = zgrid
    label1 = "Y (km)"
    label2 = "Z (km)"
if normalDirection==1:
    plotdata_2D = plotdata[:,normalIndex,:]
    grid1 = xgrid
    grid2 = zgrid
    label1 = "X (km)"
    label2 = "Z (km)"
if normalDirection==2:
    plotdata_2D = plotdata[:,:,normalIndex]
    grid1 = xgrid
    grid2 = ygrid
    label1 = "X (km)"
    label2 = "Y (km)"

# set up the plot
plt.axes().set_aspect('equal', 'datalim')
cax = plt.pcolormesh(grid1,grid2,plotdata_2D)
plt.axes().set_xlabel(label1)
plt.axes().set_ylabel(label2)
plt.colorbar(cax)
plt.show()
    
