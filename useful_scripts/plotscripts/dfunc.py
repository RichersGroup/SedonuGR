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
normalIndex = 0 # index of zone in normalDirection direction

# M0: must be 0
# M1: 0:Fx 1:Fy 2:Fz
# M2: 0:Pxx 1:Pyx 2:Pzx 3:Pyy 4:Pzy 5:Pzz
# polar: x->theta  y->phi  z->r
quantity = "M0" # "M0","M1","M2"
index = 0

##############
# NOT INPUTS #
##############

# open the hdf5 file
f = h5py.File(dataFile,"r")

# read the grid structure and erad data
xgrid = np.array(f["grid_0(cm)"])/1e5
ygrid = np.array(f["grid_1(cm)"])/1e5
zgrid = np.array(f["grid_2(cm)"])/1e5
if group>0:
    M0       = np.array(f["/distribution(erg|ccm,lab)/M0"]       )[:,:,:,species,group,0]
    plotdata = np.array(f["/distribution(erg|ccm,lab)/"+quantity])[:,:,:,species,group,index]
else:
    M0       = np.array(f["/distribution(erg|ccm,lab)/M0"]       )[:,:,:,species,:,0]
    plotdata = np.array(f["/distribution(erg|ccm,lab)/"+quantity])[:,:,:,species,:,index]    
    M0       = np.sum(M0,      axis=3)
    plotdata = np.sum(plotdata,axis=3)
    
if quantity != "M0":
    plotdata = plotdata/M0
else:
    plotdata = log10(plotdata)

print(np.min(M0),np.max(M0))
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
    
