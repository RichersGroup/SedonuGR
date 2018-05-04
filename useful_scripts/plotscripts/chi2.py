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
nfiles = 10
fileroot = "fluid_"

qlist = [
    "Pxx",
#    "Pyy",
#    "Pzz",
    "Pxy",
#    "Pxz",
#    "Pyz",
#    "Lxxx",
#    "Lxxy",
#    "Lxxz",
#    "Lxyy",
#    "Lxyz",
#    "Lxzz",
#    "Lyyy",
#    "Lyyz",
#    "Lyzz",
#    "Lzzz"
    ]
q_dof = 5
#ndof = 12 # all
#ndof = 5  # P
#ndof = 7  # L
#ndof = 2  # Pii (diag)
#ndof = 3  # Pik (off-diag)
#ndof = 3  # Liii (3diag)

#####################
# CLOSURE FUNCTIONS #
#####################
def append_data_theory(this_data, this_theory, dset, ig, q):
    E = dset[:,:,:,ig,0]
    F = dset[:,:,:,ig,1:4]
    if(q=="Pxx"):
        this_data.append(dset[:,:,:,ig,4])
        this_theory.append(E/3.)
    if(q=="Pxy"):
        this_data.append(dset[:,:,:,ig,5])
        this_theory.append(E*0.)
    if(q=="Pxz"):
        this_data.append(dset[:,:,:,ig,6])
        this_theory.append(E*0.)
    if(q=="Pyy"):
        this_data.append(dset[:,:,:,ig,7])
        this_theory.append(E/3.)
    if(q=="Pyz"):
        this_data.append(dset[:,:,:,ig,8])
        this_theory.append(E*0.)
    if(q=="Pzz"):
        this_data.append(dset[:,:,:,ig,9])
        this_theory.append(E/3.)
    if(q=="Lxxx"):
        this_data.append(dset[:,:,:,ig,10])
        this_theory.append(F[:,:,:,0]/3.)
    if(q=="Lxxy"):
        this_data.append(dset[:,:,:,ig,11])
        this_theory.append(F[:,:,:,1]/3.)
    if(q=="Lxxz"):
        this_data.append(dset[:,:,:,ig,12])
        this_theory.append(F[:,:,:,2]/3.)
    if(q=="Lxyy"):
        this_data.append(dset[:,:,:,ig,13])
        this_theory.append(F[:,:,:,0]/3.)
    if(q=="Lxyz"):
        this_data.append(dset[:,:,:,ig,14])
        this_theory.append(E[:,:,:]*0.  )
    if(q=="Lxzz"):
        this_data.append(dset[:,:,:,ig,15])
        this_theory.append(F[:,:,:,0]/3.)
    if(q=="Lyyy"):
        this_data.append(dset[:,:,:,ig,16])
        this_theory.append(F[:,:,:,1]/3.)
    if(q=="Lyyz"):
        this_data.append(dset[:,:,:,ig,17])
        this_theory.append(F[:,:,:,2]/3.)
    if(q=="Lyzz"):
        this_data.append(dset[:,:,:,ig,18])
        this_theory.append(F[:,:,:,1]/3.)
    if(q=="Lzzz"):
        this_data.append(dset[:,:,:,ig,19])
        this_theory.append(F[:,:,:,2]/3.)
        

###############
# SELECT DATA #
###############
f = h5py.File(fileroot+"00000.h5","r")
xgrid = np.array(f["axes/x0(cm)[edge]"])
ygrid = np.array(f["axes/x1(cm)[edge]"])
zgrid = np.array(f["axes/x2(cm)[edge]"])
Egrid = np.array(f["axes/frequency(Hz)[mid]"])

# neutrino energy density
for iq,q in zip(range(len(qlist)), qlist):
    print(q)
    this_dof = 0
    for fi in range(1,nfiles+1):
        this_data = []
        this_theory = []
        filename = fileroot+("%05d"%fi)+".h5"
        print(" ",filename)
        f = h5py.File(filename)

        for s in range(3):
            #print("  s="+str(s))
            dset = np.array(f["distribution"+str(s)+"(erg|ccm,tet)"])
            for ig in range(len(Egrid)):
                #print("    ig="+str(ig))
                append_data_theory(this_data, this_theory, dset, ig, q)
                this_dof += 1

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

    data     /= nfiles
    theory   /= nfiles
    data2    /= nfiles
    theory2  /= nfiles
    this_dof /= nfiles
    stddev2 = (data2 - data**2) + (theory2 - theory**2)
    stddev2[np.where(stddev2==0)] = -1.

    # report quality
    print("  this_dof: ",this_dof)
    nx = np.shape(data)
    print("  shape: ",nx)
    badGridCell = np.full(np.shape(data[0]), False, dtype=bool)
    for igs in range(nx[0]): # looping over energies and species
        badGridCell = np.logical_or(badGridCell, stddev2[igs]<0 )
    wherebad = np.where(badGridCell)
    wheregood = np.where(np.logical_not(badGridCell) )
    indices = np.array(wherebad)
    ir = sqrt(np.sum(indices*indices,axis=0))
    print("  max index radius: ",np.max(ir))
    print("  min index radius: ",np.min(ir))
    print(" ",len(ir)," bad points")

    # get chi2
    chi2 = np.sum((data-theory)**2 / stddev2, axis=0) / this_dof
    min_index = [wheregood[i][np.argmin(chi2[wheregood])] for i in range(3)]
    max_index = [wheregood[i][np.argmax(chi2[wheregood])] for i in range(3)]
    print("  min X2: ",np.min(chi2[min_index[0],min_index[1],min_index[2]]),min_index)
    print("  max X2: ",np.min(chi2[max_index[0],max_index[1],max_index[2]]),max_index)    
    avg_chi2 = np.sum(chi2[wheregood]) / np.size(chi2[wheregood])
    print("  avg X2: ",avg_chi2)

    # uncertainty-weighted relative error
    UWRE_cell = np.sum( np.abs(data-theory)/(np.abs(data)+np.abs(theory)), axis=0)
    weights_cell = np.sum(np.ones(np.shape(data)),axis=0)
    UWRE = np.sum(UWRE_cell[wheregood])
    weights = np.sum(weights_cell[wheregood])
    print("  UWRE: ",UWRE/weights)
    
    # accumulate onto net chi2
    if(iq==0):
        sum_chi2 = chi2 * this_dof
        sum_badGridCell = badGridCell
        sum_UWRE = UWRE_cell
        sum_weights = weights_cell
    else:
        sum_chi2 += chi2 * this_dof
        sum_badGridCell = np.logical_or(sum_badGridCell, badGridCell)
        sum_UWRE += UWRE_cell
        sum_weights += weights_cell
    print()

print("============================")
print("= Total over all variables =")
print("============================")

# report quality
tot_dof = this_dof * q_dof
print("q_dof: ",q_dof)
print("tot_dof: ",tot_dof)
wherebad  = np.where(sum_badGridCell)
wheregood = np.where(np.logical_not(sum_badGridCell))
indices = np.array(wherebad)
ir = sqrt(np.sum(indices*indices,axis=0))
print("max index radius: ",np.max(ir))
print("min index radius: ",np.min(ir))
print(len(ir)," bad points")

# compute total chi2
sum_chi2 /= tot_dof
min_index = [wheregood[i][np.argmin(sum_chi2[wheregood])] for i in range(3)]
max_index = [wheregood[i][np.argmax(sum_chi2[wheregood])] for i in range(3)]
print("min X2: ",np.min(sum_chi2[min_index[0],min_index[1],min_index[2]]),min_index)
print("max X2: ",np.min(sum_chi2[max_index[0],max_index[1],max_index[2]]),max_index)
avg_chi2 = np.sum(sum_chi2[wheregood]) / np.size(sum_chi2[wheregood])
print("avg X2: ",avg_chi2)

# uncertainty-weighted relative error
UWRE = np.sum(sum_UWRE[wheregood])
weights = np.sum(sum_weights[wheregood])
print("UWRE: ",UWRE/weights)

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

