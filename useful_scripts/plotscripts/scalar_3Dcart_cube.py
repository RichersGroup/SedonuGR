import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import h5py
import sys
import numpy as np
from pylab import *

##########
# INPUTS #
##########
dataFile = "fluid_combined.h5"
normalDirection = 2 # which direction is normal to the plotting plane? 0:x 1:y 2:z
normalIndex = 0 # index of zone in normalDirection direction

###############
# SELECT DATA #
###############
f = h5py.File(dataFile,"r")
xgrid = np.array(f["axes/x0(cm)[edge]"]) / 1e5
ygrid = np.array(f["axes/x1(cm)[edge]"]) / 1e5
zgrid = np.array(f["axes/x2(cm)[edge]"]) / 1e5

# gas temperature
#plotdata = log10(np.array(f["T_gas(MeV,tet)"]))

# gas density
#plotdata = np.array(f["rho(g|ccm,tet)"])

# volume factor
#plotdata = np.array(f["sqrtdetg3"])

# lapse
#plotdata = np.array(f["lapse"])

# velocity
#plotdata = np.array(f["sqrt_vdotv(cm|s,lab)"])/2.99e10

# gtt
shiftup = np.array(f["shiftup"])
metric = np.array(f["threemetric"])
g3 = np.zeros((len(xgrid)-1, \
               len(ygrid)-1, \
               len(zgrid)-1, \
               3,3))
plotdata = np.zeros((len(xgrid)-1, \
                     len(ygrid)-1, \
                     len(zgrid)-1))
data = np.array(f["threemetric"])
g3[:,:,:,0,0] = data[:,:,:,0]
g3[:,:,:,1,1] = data[:,:,:,1]
g3[:,:,:,2,2] = data[:,:,:,2]
g3[:,:,:,0,1] = g3[:,:,:,1,0] = data[:,:,:,3]
g3[:,:,:,0,2] = g3[:,:,:,2,0] = data[:,:,:,4]
g3[:,:,:,1,2] = g3[:,:,:,2,1] = data[:,:,:,5]
for i in range(0,3):
    for j in range(0,3):
        plotdata[:,:,:] = plotdata[:,:,:] + shiftup[:,:,:,i]*g3[:,:,:,i,j]*shiftup[:,:,:,j]
plotdata = plotdata - np.array(f["lapse"])

# leptonization
#data1 = np.array(f["l_emit(1|s|ccm,tet)"])
#data2 = np.array(f["l_abs(1|s|ccm,tet)"])
#plotdata = data2 - data1

# net heating
#data1 = np.array(f["four-force[abs](erg|s|g,tet)"][:,:,:,3])
#data2 = np.array(f["four-force[emit](erg|s|g,tet)"][:,:,:,3])
#plotdata = data2

# annihilation rate
#plotdata = np.array(f["annihilation_4force(erg|ccm|s,tet)"][:,:,:,3]) / np.array(f["rho(g|ccm,tet)"])

# neutrino energy density
#data0 = np.array(f["distribution0(erg|ccm,tet)"][:,:,:,:,0]).sum(axis=3)
#data1 = np.array(f["distribution1(erg|ccm,tet)"][:,:,:,:,0]).sum(axis=3)
#data2 = np.array(f["distribution2(erg|ccm,tet)"][:,:,:,:,0]).sum(axis=3)
#plotdata = data0 + data1 + data2

# neutrino shear stress
#Ptensor = np.zeros((len(xgrid)-1, \
#                    len(ygrid)-1, \
#                    len(zgrid)-1, \
#                    3,3))
#data = np.abs(np.array(f["distribution0(erg|ccm,tet)"][:,:,:,:,4:10]).sum(axis=3)) + \
#       np.abs(np.array(f["distribution1(erg|ccm,tet)"][:,:,:,:,4:10]).sum(axis=3)) + \
#       np.abs(np.array(f["distribution2(erg|ccm,tet)"][:,:,:,:,4:10]).sum(axis=3))
#Ptensor[:,:,:,0,0] = data[:,:,:,0]
#Ptensor[:,:,:,1,1] = data[:,:,:,3]
#Ptensor[:,:,:,2,2] = data[:,:,:,5]
#Ptensor[:,:,:,0,1] = Ptensor[:,:,:,1,0] = data[:,:,:,1]
#Ptensor[:,:,:,0,2] = Ptensor[:,:,:,2,0] = data[:,:,:,2]
#Ptensor[:,:,:,1,2] = Ptensor[:,:,:,2,1] = data[:,:,:,4]
#print(np.shape(Ptensor))
#eigenvals = np.array([[[ linalg.eigvalsh(Ptensor[i,j,k]) for k in range(len(zgrid)-1) ] for j in range(len(ygrid)-1) ] for i in range(len(xgrid)-1) ])
#print(np.shape(eigenvals))
#I1 = eigenvals[:,:,:,0]+eigenvals[:,:,:,1]+eigenvals[:,:,:,2]
#I2 = eigenvals[:,:,:,0]*eigenvals[:,:,:,1]+\
#     eigenvals[:,:,:,0]*eigenvals[:,:,:,2]+\
#     eigenvals[:,:,:,1]*eigenvals[:,:,:,2]
#I3 = eigenvals[:,:,:,0]*eigenvals[:,:,:,1]*eigenvals[:,:,:,2]
#plotdata = I2/(I1*I1)

# correct and take log where val is less than zero
#plotdata[np.where(I1==0)] = 0
print(np.min(plotdata),np.max(plotdata))
#minval = np.min(plotdata[np.where(plotdata>0)])
#plotdata[np.where(plotdata<=0)] = minval
#plotdata = np.log10(plotdata)
print(np.min(plotdata),np.max(plotdata))


##############
# NOT INPUTS #
##############

# make colors range from 0-1
minval = np.min(plotdata)
maxval = np.max(plotdata)
colordata = plotdata - minval
colordata = colordata / (maxval-minval)

# set up the plot
fig = plt.figure()
ax = fig.add_subplot(111,projection='3d')
ax.view_init(30,45)

cset = [[],[],[]]
Y,X = np.meshgrid(xgrid,ygrid)
Z=0*X#np.sin(np.sqrt(X*X + Y*Y)/10e5)
cset[0] = ax.plot_surface(X,Y,Z,facecolors=cm.jet(colordata[:,:,0]),rstride=1,cstride=1,linewidth=0,antialiased=False,shade=False)

Z,X = np.meshgrid(xgrid,zgrid)
Y=0*X#np.sin(np.sqrt(X*X + Y*Y)/10e5)
cset[1] = ax.plot_surface(X,Y,Z,facecolors=cm.jet(colordata[:,0,:]),rstride=1,cstride=1,linewidth=0,antialiased=False,shade=False)

Z,Y = np.meshgrid(ygrid,zgrid)
X=0*Y#np.sin(np.sqrt(X*X + Y*Y)/10e5)
cset[2] = ax.plot_surface(X,Y,Z,facecolors=cm.jet(colordata[0,:,:]),rstride=1,cstride=1,linewidth=0,antialiased=False,shade=False)

m = cm.ScalarMappable(cmap=cm.jet)
m.set_array(plotdata)
cb = fig.colorbar(m,aspect=15)
cb.set_label(r'Annihilation Rate (erg/s/g)')

ax.set_xlabel("X (km)")
ax.set_ylabel("Y (km)")
ax.set_zlabel("Z (km)")
ax.set_frame_on(False)
ax.grid(linewidth=0)

ax.plot3D([np.min(xgrid),np.max(xgrid)],[0,0],[0,0],color='k',linewidth=1.5)
ax.plot3D([0,0],[np.min(ygrid),np.max(ygrid)],[0,0],color='k',linewidth=1.5)
ax.plot3D([0,0],[0,0],[np.min(zgrid),np.max(zgrid)],color='k',linewidth=1.5)

# set up colorbar
#plt.colorbar()

#plt.show()

plt.savefig("threeslice.png",bbox_inches="tight",dpi=150)
