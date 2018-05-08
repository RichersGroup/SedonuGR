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
nfiles = 10
interpolator = "thick"
species_list = [0,1,2]
igmin = 0
igmax = 16

qlist = [
    "Pxx",
    "Pyy",
    "Pzz",
    "Pxy",
    "Pxz",
    "Pyz",
    "Lxxx",
    "Lxxy",
    "Lxxz",
    "Lxyy",
    "Lxyz",
    "Lxzz",
    "Lyyy",
    "Lyyz",
    "Lyzz",
    "Lzzz"
    ]
q_dof = 12
#ndof = 12 # all
#ndof = 5  # P
#ndof = 7  # L
#ndof = 2  # Pii (diag)
#ndof = 3  # Pik (off-diag)
#ndof = 3  # Liii (3diag)

qdict = {
    "E" : 0,
    "Fx": 1,
    "Fy": 2,
    "Fz": 3,
    "Pxx": 4,
    "Pxy": 5,
    "Pxz": 6,
    "Pyy": 7,
    "Pyz": 8,
    "Pzz": 9,
    "Lxxx": 10,
    "Lxxy": 11,
    "Lxxz": 12,
    "Lxyy": 13,
    "Lxyz": 14,
    "Lxzz": 15,
    "Lyyy": 16,
    "Lyyz": 17,
    "Lyzz": 18,
    "Lzzz": 19
}

#####################
# CLOSURE FUNCTIONS #
#####################
def q_thick_thin(E, Fx, Fy, Fz, FdotF, q):
    if(q=="Pxx"):
        return E/3.                 , E * Fx*Fx/FdotF
    if(q=="Pxy"):
        return np.zeros(np.shape(E)), E * Fx*Fy/FdotF
    if(q=="Pxz"):
        return np.zeros(np.shape(E)), E * Fx*Fz/FdotF
    if(q=="Pyy"):
        return E/3.                 , E * Fy*Fy/FdotF
    if(q=="Pyz"):
        return np.zeros(np.shape(E)), E * Fy*Fz/FdotF
    if(q=="Pzz"):
        return E/3.                 , E * Fz*Fz/FdotF
    if(q=="Lxxx"):
        return Fx/3.                , Fx * Fx*Fx/FdotF
    if(q=="Lxxy"):
        return Fy/3.                , Fy * Fx*Fx/FdotF
    if(q=="Lxxz"):
        return Fz/3.                , Fz * Fx*Fx/FdotF
    if(q=="Lxyy"):
        return Fx/3.                , Fx * Fy*Fy/FdotF
    if(q=="Lxyz"):
        return np.zeros(np.shape(E)), np.zeros(np.shape(E))
    if(q=="Lxzz"):
        return Fx/3.                , Fx * Fz*Fz/FdotF
    if(q=="Lyyy"):
        return Fy/3.                , Fy * Fy*Fy/FdotF
    if(q=="Lyyz"):
        return Fz/3.                , Fz * Fy*Fy/FdotF
    if(q=="Lyzz"):
        return Fy/3.                , Fy * Fz*Fz/FdotF
    if(q=="Lzzz"):
        return Fz/3.                , Fz * Fz*Fz/FdotF

def closure(dset, ig, q):
    E  = dset[:,:,:,ig,qdict["E"]]
    Fx = dset[:,:,:,ig,qdict["Fx"]]
    Fy = dset[:,:,:,ig,qdict["Fy"]]
    Fz = dset[:,:,:,ig,qdict["Fz"]]
    FdotF = Fx*Fx + Fy*Fy + Fz*Fz
    FdotF[np.where(FdotF==0)] = math.inf # force free-streaming fluxfac to be zero
    fluxfac = np.sqrt(FdotF) / E
    fluxfac[np.where(FdotF==0)] = 0.0
    #print("fluxfac max=",np.max(fluxfac)) # should be 1 at most
    #print("fluxfac min=",np.min(fluxfac)) # should be 0 at least

    if interpolator=="thick":
        chi = np.ones(np.shape(E)) / 3.
    if interpolator=="thin":
        chi = np.ones(np.shape(E))
    if interpolator=="Kershaw":
        chi = (1. + 2.*fluxfac) / 3.
    if interpolator=="Wilson":
        chi = 1./3. - 1./3.*fluxfac + fluxfac**2
    if interpolator=="Levermore":
        chi = (3 + 4.*fluxfac**2) / (5. + 2.*np.sqrt(4. - 3.*fluxfac**2))
    if interpolator=="MEFD_maxpack":
        chi = 1./3. * (1.0 -2.*fluxfac + 4.*fluxfac**2)
    if interpolator=="ME":
        chi = 1./3. + (2.*fluxfac**2)/15. * (3. - fluxfac + 3.*fluxfac**2)
    if interpolator=="Janka1":
        a=0.5
        b=1.3064
        n=4.1342
        chi = 1./3. * (1. + a*fluxfac**m + (2.-a)*fluxfac**n)
    if interpolator=="Janka2":
        a=1.0
        b=1.345
        n=5.1717
        chi = 1./3. * (1. + a*fluxfac**m + (2.-a)*fluxfac**n)
    #print("chi max=",np.max(chi)) # should be 1 at most
    #print("chi min=",np.min(chi)) # should be 1/3 at least

    q_thick, q_thin = q_thick_thin(E, Fx, Fy, Fz, FdotF, q)
    return ( 3.*chi-1.0  )/2. * q_thin \
        +  ( 3.*(1.-chi) )/2. * q_thick
    


###############
# SELECT DATA #
###############
fileroot = "fluid_"
f = h5py.File(fileroot+"00000.h5","r")
xgrid = np.array(f["axes/x0(cm)[edge]"])
ygrid = np.array(f["axes/x1(cm)[edge]"])
zgrid = np.array(f["axes/x2(cm)[edge]"])
Egrid = np.array(f["axes/frequency(Hz)[mid]"])
nx = len(xgrid)-1
ny = len(ygrid)-1
nz = len(zgrid)-1
ng = len(Egrid[igmin:igmax])
nvars_per_q = ng*len(species_list)

# neutrino energy density
print("igmin =",igmin)
print("igmax =",igmax)
print("species_list =",species_list)
print("nfiles =",nfiles)
print("interpolator =",interpolator)
print()

allq_X2 = np.zeros((nx,ny,nz))
allq_A2 = np.zeros((nx,ny,nz))
allq_badGridCell = np.full((nx,ny,nz), False, dtype=bool)
for iq,q in zip(range(len(qlist)), qlist):
    print(q)
    data    = np.zeros((nvars_per_q,nx,ny,nz))
    data2   = np.zeros((nvars_per_q,nx,ny,nz))
    theory  = np.zeros((nvars_per_q,nx,ny,nz))
    theory2 = np.zeros((nvars_per_q,nx,ny,nz))
    for fi in range(1,nfiles+1):
        filename = fileroot+("%05d"%fi)+".h5"
        #print(" ",filename)
        f = h5py.File(filename)
        this_data   = np.zeros((nvars_per_q,nx,ny,nz))
        this_theory = np.zeros((nvars_per_q,nx,ny,nz))

        for s in species_list:
            dset = np.array(f["distribution"+str(s)+"(erg|ccm,tet)"])
            for ig in range(igmin, igmax):
                index = ig + s*ng
                this_data[index,:,:,:]   = dset[:,:,:,ig,qdict[q]]
                this_theory[index,:,:,:] = closure(dset, ig, q)

        data    += this_data
        data2   += this_data**2
        theory  += this_theory
        theory2 += this_theory**2
        f.close()

    data     /= nfiles
    theory   /= nfiles
    data2    /= nfiles
    theory2  /= nfiles
    stddev2 = ( (data2 - data**2) + (theory2 - theory**2) ) / nfiles # stddev of mean
    stddev2[np.where(stddev2==0)] = -1.

    # report quality
    print("  shape:",np.shape(data))
    badGridCell = np.full((nx,ny,nz), False, dtype=bool)
    for igs in range(nvars_per_q): # looping over energies and species
        badGridCell = np.logical_or(badGridCell, stddev2[igs]<0 )
    wherebad = np.where(badGridCell)
    wheregood = np.where(np.logical_not(badGridCell) )
    indices = np.array(wherebad)
    ir = sqrt(np.sum(indices*indices,axis=0))
    print("  max index radius:",np.max(ir))
    print("  min index radius:",np.min(ir))
    print(" ",len(ir)," bad points")
    ngood = np.shape(wheregood)[1]
    
    # sum over species and energy bins
    X2 = np.sum((       data -       theory )**2 / stddev2 , axis=0)
    A2 = np.sum((np.abs(data)+np.abs(theory))**2 / stddev2 , axis=0)
    min_index = [wheregood[i][np.argmin(X2[wheregood])] for i in range(3)]
    max_index = [wheregood[i][np.argmax(X2[wheregood])] for i in range(3)]
    print("  min cell-summed X2/dof:",X2[min_index[0],min_index[1],min_index[2]]/nvars_per_q, min_index)
    print("  max cell-summed X2/dof:",X2[max_index[0],max_index[1],max_index[2]]/nvars_per_q, max_index)
    min_index = [wheregood[i][np.argmin(A2[wheregood])] for i in range(3)]
    max_index = [wheregood[i][np.argmax(A2[wheregood])] for i in range(3)]
    print("  min cell-summed A2/dof:",A2[min_index[0],min_index[1],min_index[2]]/nvars_per_q, min_index)
    print("  max cell-summed A2/dof:",A2[max_index[0],max_index[1],max_index[2]]/nvars_per_q, max_index)

    # sum over cell to get scalar
    X2_dof = np.sum(X2[wheregood]) / (ngood*nvars_per_q)
    A2_dof = np.sum(A2[wheregood]) / (ngood*nvars_per_q)
    print("  total X2/dof:", X2_dof)
    print("  total A2/dof:", A2_dof)
    print("  avg rel_error:", sqrt(X2_dof/A2_dof))

    # accumulate onto net chi2
    allq_X2 += X2
    allq_A2 += A2
    allq_badGridCell = np.logical_or(allq_badGridCell, badGridCell)
    print()

print("============================")
print("= Total over all variables =")
print("============================")

# report quality
tot_dof = nvars_per_q * q_dof
print("q_dof:",q_dof)
print("tot_dof:",tot_dof)
wherebad  = np.where(allq_badGridCell)
wheregood = np.where(np.logical_not(allq_badGridCell))
indices = np.array(wherebad)
ir = sqrt(np.sum(indices*indices,axis=0))
print("max index radius:",np.max(ir))
print("min index radius:",np.min(ir))
print(len(ir)," bad points")
ngood = np.shape(wheregood)[1]

# compute total chi2
min_index = [wheregood[i][np.argmin(allq_X2[wheregood])] for i in range(3)]
max_index = [wheregood[i][np.argmax(allq_X2[wheregood])] for i in range(3)]
print("min X2/dof:",allq_X2[min_index[0],min_index[1],min_index[2]]/tot_dof, min_index)
print("max X2/dof:",allq_X2[max_index[0],max_index[1],max_index[2]]/tot_dof, max_index)
min_index = [wheregood[i][np.argmin(allq_A2[wheregood])] for i in range(3)]
max_index = [wheregood[i][np.argmax(allq_A2[wheregood])] for i in range(3)]
print("min A2/dof:",allq_A2[min_index[0],min_index[1],min_index[2]]/tot_dof, min_index)
print("max A2/dof:",allq_A2[max_index[0],max_index[1],max_index[2]]/tot_dof, max_index)

# sum over cell to get scalar
X2_dof = np.sum(allq_X2[wheregood]) / (ngood*tot_dof)
A2_dof = np.sum(allq_A2[wheregood]) / (ngood*tot_dof)
print("total X2/dof:", X2_dof)
print("total A2/dof:", A2_dof)
print("avg rel_error:", sqrt(X2_dof/A2_dof))

exit()

##############
# NOT INPUTS #
##############
normalDirection = 0 # which direction is normal to the plotting plane? 0:x 1:y 2:z
normalIndex = 0 # index of zone in normalDirection direction

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

