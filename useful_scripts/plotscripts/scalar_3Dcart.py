import matplotlib.pyplot as plt
import h5py
import sys
import numpy as np
from pylab import *

##########
# INPUTS #
##########
dataFile = "/mnt/scratch/Christian_M1_Simulations/3D/fluid_00000.h5"
species = 0 # 0:nue 1:nuebar 2:nux
normalDirection = 2 # which direction is normal to the plotting plane? 0:x 1:y 2:z
normalIndex = 0 # index of zone in normalDirection direction

##############
# NOT INPUTS #
##############

# open the hdf5 file
f = h5py.File(dataFile,"r")

# read the grid structure and erad data
xgrid = np.array(f["grid_0(cm)"])
ygrid = np.array(f["grid_1(cm)"])
zgrid = np.array(f["grid_2(cm)"])
#plotdata = log10(np.array(f["e_rad(erg|ccm,lab)"])[:,:,:,species])
plotdata = log10(np.array(f["T_gas(MeV,com)"]))
#plotdata = log10(np.array(f["C_emit(erg|s|g,com)"]))
#plotdata = log10(np.array(f["H_abs(erg|s|g,com)"]))
#plotdata = log10(np.abs(np.array(f["dYe_dt_abs(1|s,lab)"])))

# select out only the plane we want
if normalDirection==0:
    plotdata_2D = plotdata[normalIndex,:,:]
    grid1 = ygrid
    grid2 = zgrid
if normalDirection==1:
    plotdata_2D = plotdata[:,normalIndex,:]
    grid1 = xgrid
    grid2 = zgrid
if normalDirection==2:
    plotdata_2D = plotdata[:,:,normalIndex]
    grid1 = xgrid
    grid2 = ygrid

# set up the plot
plt.axes().set_aspect('equal', 'datalim')
plt.pcolormesh(grid1,grid2,plotdata_2D)
plt.show()
    
