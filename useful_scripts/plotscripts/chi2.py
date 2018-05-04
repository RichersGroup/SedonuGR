import matplotlib.pyplot as plt
import h5py
import sys
import numpy as np
from pylab import *
import random
random.seed()
##########
# INPUTS #
##########
normalDirection = 0 # which direction is normal to the plotting plane? 0:x 1:y 2:z
normalIndex = 0 # index of zone in normalDirection direction
nfiles = 3
fileroot = "fluid_"

#####################
# CLOSURE FUNCTIONS #
#####################
dof_P = 5
def append_P_data(this_data, f, dset, ig):
    this_data.append(dset[:,:,:,ig,4]) # xx
    this_data.append(dset[:,:,:,ig,5]) # xy
    this_data.append(dset[:,:,:,ig,6]) # xz
    this_data.append(dset[:,:,:,ig,7]) # yy
    this_data.append(dset[:,:,:,ig,8]) # yz
    this_data.append(dset[:,:,:,ig,9]) # zz
def append_P_theory(this_theory, f, dset, ig):
    E = dset[:,:,:,ig,0]
    #F = np.array(f[dsetname][:,:,:,ig,1:4])
    this_theory.append(E/3.) # xx
    this_theory.append(E*0.) # xy
    this_theory.append(E*0.) # xz
    this_theory.append(E/3.) # yy
    this_theory.append(E*0.) # yz
    this_theory.append(E/3.) # zz

dof_L = 7
def append_L_data(this_data, f, dset, ig):
    this_data.append(dset[:,:,:,ig,10]) # xxx
    this_data.append(dset[:,:,:,ig,11]) # xxy
    this_data.append(dset[:,:,:,ig,12]) # xxz
    this_data.append(dset[:,:,:,ig,13]) # xyy
    this_data.append(dset[:,:,:,ig,14]) # xyz
    this_data.append(dset[:,:,:,ig,15]) # xzz
    this_data.append(dset[:,:,:,ig,16]) # yyy
    this_data.append(dset[:,:,:,ig,17]) # yyz
    this_data.append(dset[:,:,:,ig,18]) # yzz
    this_data.append(dset[:,:,:,ig,19]) # zzz
def append_L_theory(this_theory, f, dset, ig):
    E = dset[:,:,:,ig,0]
    F = dset[:,:,:,ig,1:4]
    this_theory.append(F[:,:,:,0]/3.) # xxx
    this_theory.append(F[:,:,:,1]/3.) # xxy
    this_theory.append(F[:,:,:,2]/3.) # xxz
    this_theory.append(F[:,:,:,0]/3.) # xyy
    this_theory.append(E[:,:,:]*0.  ) # xyz
    this_theory.append(F[:,:,:,0]/3.) # xzz
    this_theory.append(F[:,:,:,1]/3.) # yyy
    this_theory.append(F[:,:,:,2]/3.) # yyz
    this_theory.append(F[:,:,:,1]/3.) # yzz
    this_theory.append(F[:,:,:,2]/3.) # zzz

###############
# SELECT DATA #
###############
f = h5py.File(fileroot+"00000.h5","r")
xgrid = np.array(f["axes/x0(cm)[edge]"])
ygrid = np.array(f["axes/x1(cm)[edge]"])
zgrid = np.array(f["axes/x2(cm)[edge]"])
Egrid = np.array(f["axes/frequency(Hz)[mid]"])

# get the metric information
#g3 = np.array(f["threemetric"])
#lapse = np.array(f["lapse"])
#tmp = np.array(f["shiftup"])
#shiftup = np.array([
#    tmp[:,:,:,0],
#    tmp[:,:,:,1],
#    tmp[:,:,:,2]
#    ])
#print(np.shape(shiftup))
#shiftdown = np.array([
#    tmp[:,:,:,0]*g3[:,:,:,0] + tmp[:,:,:,1]*g3[:,:,:,3] + tmp[:,:,:,2]*g3[:,:,:,4],
#    tmp[:,:,:,0]*g3[:,:,:,3] + tmp[:,:,:,1]*g3[:,:,:,1] + tmp[:,:,:,2]*g3[:,:,:,5],
#    tmp[:,:,:,0]*g3[:,:,:,4] + tmp[:,:,:,1]*g3[:,:,:,5] + tmp[:,:,:,2]*g3[:,:,:,2]
#    ])
#print(np.shape(shiftdown))
#gtt = -lapse + sum(shiftup*shiftdown,axis=0)
#wherebad = np.where(gtt>0)

# neutrino energy density
for fi in range(1,nfiles+1):
    this_data = []
    this_theory = []
    filename = fileroot+("%05d"%fi)+".h5"
    print(filename)
    f = h5py.File(filename)

    for s in range(3):
        #print("  s="+str(s))
        dset = np.array(f["distribution"+str(s)+"(erg|ccm,tet)"])
        for ig in range(len(Egrid)):
            #print("    ig="+str(ig))
            
            # neutrino pressure tensor
            #dof += dof_P
            append_P_data(    this_data, f, dset, ig)
            append_P_theory(this_theory, f, dset, ig)

            # neutrino L tensor
            #dof += dof_L
            append_L_data(    this_data, f, dset, ig)
            append_L_theory(this_theory, f, dset, ig)
            
    this_data = np.array(this_data)
    this_theory = np.array(this_theory)
    if fi==1:
        data    = this_data
        data2   = this_data**2
        theory  = this_theory
        theory2 = this_theory**2
    else:
        data    += this_data
        data2   += this_data**2
        theory  += this_theory
        theory2 += this_theory**2
    f.close()

data    /= nfiles
theory  /= nfiles
data2   /= nfiles
theory2 /= nfiles
data_stddev2   =   data2 -   data**2
theory_stddev2 = theory2 - theory**2
diff2_stddev2 = (data-theory)**2 / (data_stddev2 + theory_stddev2)
diff2_stddev2[np.where(diff2_stddev2!=diff2_stddev2)] = 0

# get the number of degrees of freedom
ndof = dof_P + dof_L #len(data[:,0,0,0])
print("ndof: ",ndof)
nx = np.shape(data[0,:,:,:])
print("shape: ",nx)
dof = np.array( [[[np.count_nonzero(diff2_stddev2[:,k,j,i]) for i in range(nx[0])] for j in range(nx[1])] for k in range(nx[2])] )
wherebad = np.where(dof < ndof)
wheregood = np.where(dof==ndof)

# get chi2
chi2 = np.sum(diff2_stddev2, axis=0) / ndof
print("min X2: ",np.min(chi2[wheregood]))
print("max X2: ",np.max(chi2), np.unravel_index(np.argmax(chi2),np.shape(chi2)))

avg_chi2 = np.sum(chi2[wheregood]) / np.size(chi2[wheregood])
print("avg X2: ",avg_chi2)

# report quality
indices = np.array(wherebad)
ir = sqrt(np.sum(indices*indices,axis=0))
print("max index radius: ",np.max(ir))
print("min index radius: ",np.min(ir))
print(len(ir)," bad points")
exit()


##############
# NOT INPUTS #
##############

# select out only the plane we want
if normalDirection==0:
    plotdata_2D = data[normalIndex,:,:]
    grid1 = zgrid
    grid2 = ygrid
if normalDirection==1:
    plotdata_2D = data[:,normalIndex,:]
    grid1 = zgrid
    grid2 = xgrid
if normalDirection==2:
    plotdata_2D = data[:,:,normalIndex]
    grid1 = ygrid
    grid2 = xgrid

plotdata_2Da = np.array([[random.gauss(plotdata_2D[j,i],plotdata_2D[j,i]/5) for i in range(np.shape(plotdata_2D)[0])] for j in range(np.shape(plotdata_2D)[1])])
plotdata_2Db = np.array([[random.gauss(plotdata_2D[j,i],plotdata_2D[j,i]/5) for i in range(np.shape(plotdata_2D)[0])] for j in range(np.shape(plotdata_2D)[1])])
plotdata_2Dc = np.array([[random.gauss(plotdata_2D[j,i],plotdata_2D[j,i]/5) for i in range(np.shape(plotdata_2D)[0])] for j in range(np.shape(plotdata_2D)[1])])

plotdata_mean  = (plotdata_2Da    + plotdata_2Db    + plotdata_2Dc   ) / 3.
plotdata_mean2 = (plotdata_2Da**2 + plotdata_2Db**2 + plotdata_2Dc**2) / 3.
plotdata_stddev = np.sqrt(plotdata_mean2 - plotdata_mean**2)

plotdata_chi2 = (plotdata_2Da - plotdata_mean)**2 + (plotdata_2Db - plotdata_mean)**2 + (plotdata_2Dc - plotdata_mean)**2
plotdata_chi2 = plotdata_chi2 / (3*plotdata_mean**2)
plotdata_chi2[np.where(plotdata_chi2 != plotdata_chi2)] = 0

asdf = plotdata_chi2[np.where(plotdata_chi2>0)]
print(np.size(asdf))
int_chi2 = np.sum(asdf) / np.size(asdf)
print(int_chi2)

# set up the plot
plt.axes().set_aspect('equal', 'datalim')
plt.pcolormesh(grid1,grid2,plotdata_chi2)

# set up colorbar
plt.colorbar()

plt.show()

