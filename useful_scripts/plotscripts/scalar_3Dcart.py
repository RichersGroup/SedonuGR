import matplotlib.pyplot as plt
import h5py
import sys
import numpy as np
from pylab import *

##########
# INPUTS #
##########
dataFile = "fluid_combined.h5"
normalDirection = 0 # which direction is normal to the plotting plane? 0:x 1:y 2:z
normalIndex = 89 # index of zone in normalDirection direction

###############
# SELECT DATA #
###############
f = h5py.File(dataFile,"r")
xgrid = np.array(f["axes/x0(cm)[edge]"])
ygrid = np.array(f["axes/x1(cm)[edge]"])
zgrid = np.array(f["axes/x2(cm)[edge]"])

# gas temperature
#plotdata = log10(np.array(f["T_gas(MeV,tet)"]))

# gas density
#plotdata = log10(np.array(f["rho(g|ccm,tet)"]))

# volume factor
#plotdata = np.array(f["sqrtdetg3"])

# lapse
#plotdata = np.array(f["lapse"])

# velocity
#plotdata = np.array(f["sqrt_vdotv(cm|s,lab)"])

# leptonization
#data1 = np.array(f["l_emit(1|s|ccm,tet)"])
#data2 = np.array(f["l_abs(1|s|ccm,tet)"])
#plotdata = data2 - data1

# net heating
#data1 = np.array(f["four-force[abs](erg|s|g,tet)"][:,:,:,3])
#data2 = np.array(f["four-force[emit](erg|s|g,tet)"][:,:,:,3])
#plotdata = data2

# annihilation rate
#plotdata = np.array(f["annihilation_4force(erg|ccm|s,tet)"][:,:,:,3])

# neutrino energy density
#data0 = np.array(f["distribution0(erg|ccm,tet)"][:,:,:,:,0]).sum(axis=3)
#data1 = np.array(f["distribution1(erg|ccm,tet)"][:,:,:,:,0]).sum(axis=3)
#data2 = np.array(f["distribution2(erg|ccm,tet)"][:,:,:,:,0]).sum(axis=3)
#plotdata = data0 + data1 + data2

# neutrino shear stress
Ptensor = np.zeros((len(xgrid)-1, \
                    len(ygrid)-1, \
                    len(zgrid)-1, \
                    3,3))
data = np.abs(np.array(f["distribution0(erg|ccm,tet)"][:,:,:,:,4:10]).sum(axis=3)) + \
       np.abs(np.array(f["distribution1(erg|ccm,tet)"][:,:,:,:,4:10]).sum(axis=3)) + \
       np.abs(np.array(f["distribution2(erg|ccm,tet)"][:,:,:,:,4:10]).sum(axis=3))
Ptensor[:,:,:,0,0] = data[:,:,:,0]
Ptensor[:,:,:,1,1] = data[:,:,:,3]
Ptensor[:,:,:,2,2] = data[:,:,:,5]
Ptensor[:,:,:,0,1] = Ptensor[:,:,:,1,0] = data[:,:,:,1]
Ptensor[:,:,:,0,2] = Ptensor[:,:,:,2,0] = data[:,:,:,2]
Ptensor[:,:,:,1,2] = Ptensor[:,:,:,2,1] = data[:,:,:,4]
print(np.shape(Ptensor))
eigenvals = np.array([[[ linalg.eigvalsh(Ptensor[i,j,k]) for k in range(len(zgrid)-1) ] for j in range(len(ygrid)-1) ] for i in range(len(xgrid)-1) ])
print(np.shape(eigenvals))
I1 = eigenvals[:,:,:,0]+eigenvals[:,:,:,1]+eigenvals[:,:,:,2]
I2 = eigenvals[:,:,:,0]*eigenvals[:,:,:,1]+\
     eigenvals[:,:,:,0]*eigenvals[:,:,:,2]+\
     eigenvals[:,:,:,1]*eigenvals[:,:,:,2]
I3 = eigenvals[:,:,:,0]*eigenvals[:,:,:,1]*eigenvals[:,:,:,2]
plotdata = I2/(I1*I1)

# correct and take log where val is less than zero
plotdata[np.where(I1==0)] = 0
minval = np.min(plotdata[np.where(plotdata>0)])
plotdata[np.where(plotdata<=0)] = minval
#plotdata = np.log10(plotdata)


##############
# NOT INPUTS #
##############

# read the grid structure and erad data

# select out only the plane we want
if normalDirection==0:
    plotdata_2D = plotdata[normalIndex,:,:]
    grid1 = zgrid
    grid2 = ygrid
if normalDirection==1:
    plotdata_2D = plotdata[:,normalIndex,:]
    grid1 = zgrid
    grid2 = xgrid
if normalDirection==2:
    plotdata_2D = plotdata[:,:,normalIndex]
    grid1 = ygrid
    grid2 = xgrid

# set up the plot
plt.axes().set_aspect('equal', 'datalim')
plt.pcolormesh(grid1,grid2,plotdata_2D)

# set up colorbar
plt.colorbar()

plt.show()

