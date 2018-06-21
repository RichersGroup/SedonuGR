import matplotlib.pyplot as plt
import h5py
import sys
import numpy as np
from pylab import *

##########
# INPUTS #
##########
nfiles = sys.argv[1]
species_list = [0,1,2]
igmin = 0
igmax = 16
nFluxFacBins = 20
interpolator = sys.argv[2]

qlist = []
for iq in range(3,len(sys.argv)):
    qlist.append(sys.argv[iq])
    
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

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

#####################
# CLOSURE FUNCTIONS #
#####################
def q_thick_thin(E, F, FdotF, q):
    
    if(q=="Pxx"):
        return E/3.                 , E * F[:,:,:,:,0]*F[:,:,:,:,0]/FdotF
    if(q=="Pxy"):
        return np.zeros(np.shape(E)), E * F[:,:,:,:,0]*F[:,:,:,:,1]/FdotF
    if(q=="Pxz"):
        return np.zeros(np.shape(E)), E * F[:,:,:,:,0]*F[:,:,:,:,2]/FdotF
    if(q=="Pyy"):
        return E/3.                 , E * F[:,:,:,:,1]*F[:,:,:,:,1]/FdotF
    if(q=="Pyz"):
        return np.zeros(np.shape(E)), E * F[:,:,:,:,1]*F[:,:,:,:,2]/FdotF
    if(q=="Pzz"):
        return E/3.                 , E * F[:,:,:,:,2]*F[:,:,:,:,2]/FdotF
    if(q=="Lxxx"):
        return F[:,:,:,:,0]*3./5.   , F[:,:,:,:,0] * F[:,:,:,:,0]*F[:,:,:,:,0]/FdotF
    if(q=="Lxxy"):
        return F[:,:,:,:,1]*1./5.   , F[:,:,:,:,1] * F[:,:,:,:,0]*F[:,:,:,:,0]/FdotF
    if(q=="Lxxz"):
        return F[:,:,:,:,2]*1./5.   , F[:,:,:,:,2] * F[:,:,:,:,0]*F[:,:,:,:,0]/FdotF
    if(q=="Lxyy"):
        return F[:,:,:,:,0]*1./5.   , F[:,:,:,:,0] * F[:,:,:,:,1]*F[:,:,:,:,1]/FdotF
    if(q=="Lxyz"):
        return np.zeros(np.shape(E)), F[:,:,:,:,0] * F[:,:,:,:,1]*F[:,:,:,:,2]/FdotF
    if(q=="Lxzz"):
        return F[:,:,:,:,0]*1./5.   , F[:,:,:,:,0] * F[:,:,:,:,2]*F[:,:,:,:,2]/FdotF
    if(q=="Lyyy"):
        return F[:,:,:,:,1]*3./5.   , F[:,:,:,:,1] * F[:,:,:,:,1]*F[:,:,:,:,1]/FdotF
    if(q=="Lyyz"):
        return F[:,:,:,:,2]*1./5.   , F[:,:,:,:,2] * F[:,:,:,:,1]*F[:,:,:,:,1]/FdotF
    if(q=="Lyzz"):
        return F[:,:,:,:,1]*1./5.   , F[:,:,:,:,1] * F[:,:,:,:,2]*F[:,:,:,:,2]/FdotF
    if(q=="Lzzz"):
        return F[:,:,:,:,2]*3./5.   , F[:,:,:,:,2] * F[:,:,:,:,2]*F[:,:,:,:,2]/FdotF

def closure(E, F, q):
    with np.errstate(invalid='ignore'):
        FdotF = np.sum(F**2, axis=4)
        with np.errstate(invalid='ignore'):
            fluxfac = np.sqrt(FdotF) / E
        #print("fluxfac max=",np.max(fluxfac[np.where(fluxfac>=0)]))
        #print("fluxfac min=",np.min(fluxfac[np.where(fluxfac>=0)]))
        
        if is_number(interpolator):
            chi = np.ones(np.shape(E)) * float(interpolator)
        elif interpolator=="Kershaw":
            chi = (1. + 2.*fluxfac**2) / 3.
        elif interpolator=="Wilson":
            chi = 1./3. - 1./3.*fluxfac + fluxfac**2
        elif interpolator=="Levermore":
            chi = (3 + 4.*fluxfac**2) / (5. + 2.*np.sqrt(4. - 3.*fluxfac**2))
        elif interpolator=="MEFD_maxpack":
            chi = 1./3. * (1.0 -2.*fluxfac + 4.*fluxfac**2)
        elif interpolator=="ME":
            chi = 1./3. + (2.*fluxfac**2)/15. * (3. - fluxfac + 3.*fluxfac**2)
        elif interpolator=="Janka1":
            a=0.5
            b=1.3064
            n=4.1342
            chi = 1./3. * (1. + a*fluxfac**m + (2.-a)*fluxfac**n)
        elif interpolator=="Janka2":
            a=1.0
            b=1.345
            n=5.1717
            chi = 1./3. * (1. + a*fluxfac**m + (2.-a)*fluxfac**n)
        else:
            print("Error - unknown interpolator")
            sys.exit()
        #print("chi max=",np.max(chi[np.where(fluxfac>=0)]))
        #print("chi min=",np.min(chi[np.where(fluxfac>=0)]))

        q_thick, q_thin = q_thick_thin(E, F, FdotF, q)
        return ( 3.*chi-1.0  )/2. * q_thin \
            +  ( 3.*(1.-chi) )/2. * q_thick
    


###############
# SELECT DATA #
###############
fileroot = "fluid_"
f = h5py.File(fileroot+"00001.h5","r")
xgrid = np.array(f["axes/x0(cm)[edge]"])
ygrid = np.array(f["axes/x1(cm)[edge]"])
zgrid = np.array(f["axes/x2(cm)[edge]"])
Egrid = np.array(f["axes/frequency(Hz)[mid]"])
nx = len(xgrid)-1
ny = len(ygrid)-1
nz = len(zgrid)-1
ng = len(Egrid[igmin:igmax])
ns = len(species_list)
nvars_per_q = ng*ns

# get output filenames
outfilename = "analysis/"
if(is_number(interpolator)):
    outfilename = outfilename + "%.3f_"%float(interpolator)
else:
    outfilename = outfilename + interpolator + "_"
outfilename = outfilename + "s"
for s in range(ns):
    outfilename = outfilename + str(s)
outfilename = outfilename + "_g"+str(igmin)+"-"+str(igmax)+"_"
outfilename = outfilename + "nf"+str(nfiles)+"_"
for iq in range(len(qlist)):
    outfilename = outfilename + qlist[iq]
ftxt = open(outfilename+".txt","w")
    
# neutrino energy density
ftxt.write("igmin = "+str(igmin)+"\n")
ftxt.write("igmax = "+str(igmax)+"\n")
ftxt.write("species_list = "+str(species_list)+"\n")
ftxt.write("nfiles = "+str(nfiles)+"\n")
ftxt.write("interpolator = "+str(interpolator)+"\n")
ftxt.write("qlist = "+ str(qlist)+"\n")
ftxt.write("\n")

allq_X2         = np.zeros((ns,nx,ny,nz,ng))
allq_delta      = np.zeros((ns,nx,ny,nz,ng))
allq_badGridCell = np.full((ns,nx,ny,nz,ng), False, dtype=bool)
Etmp = np.zeros((ns,nx,ny,nz,ng))
Ftmp = np.zeros((ns,nx,ny,nz,ng,3))
for iq,q in zip(range(len(qlist)), qlist):
    ftxt.write(q+"\n")
    data    = np.zeros((ns,nx,ny,nz,ng))
    data2   = np.zeros((ns,nx,ny,nz,ng))
    theory  = np.zeros((ns,nx,ny,nz,ng))
    theory2 = np.zeros((ns,nx,ny,nz,ng))
    for fi in range(1,nfiles+1):
        filename = fileroot+("%05d"%fi)+".h5"
        ftxt.write("  "+filename+"\n")
        ftxt.flush()
        f = h5py.File(filename)
        this_data   = np.zeros((ns,nx,ny,nz,ng))
        this_theory = np.zeros((ns,nx,ny,nz,ng))

        for s in species_list:
            dsetname = "distribution"+str(s)+"(erg|ccm,tet)"
            this_data[s,:,:,:,:] = np.array(f[dsetname][:,:,:,igmin:igmax,qdict[q]])
            E = np.array(f[dsetname][:,:,:,igmin:igmax,0]) 
            F = np.array(f[dsetname][:,:,:,igmin:igmax,1:4]) 
            this_theory[s,:,:,:,:] = closure(E, F, q)
            if(iq==0):
                Etmp[s] += E
                Ftmp[s] += F
            
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
    if(iq==0):
        Etmp /= nfiles
        Ftmp /= nfiles
        with np.errstate(invalid='ignore'):
            fluxfac = np.sqrt(np.sum(Ftmp*Ftmp,axis=5)) / Etmp

    
    # report quality
    with np.errstate(invalid='ignore'):
        badGridCell = np.logical_not(stddev2>0) # catches nans too
    wherebad = np.where(badGridCell)
    wheregood = np.where(np.logical_not(badGridCell))
    ngood = len(wheregood[0])
    ftxt.write("  "+str(len(wherebad[0]))+" bad points"+"\n")
    
    # get chi2 and delta for this variable
    with np.errstate(invalid='ignore'):
        X2    = (data-theory)**2 / stddev2
        delta = np.abs(data-theory) / Etmp#/ (np.abs(data)+np.abs(theory))
    ftxt.write("  total X2/dof: "   +str(np.sum(   X2[wheregood]) / ngood)+"\n")
    ftxt.write("  total delta/dof: "+str(np.sum(delta[wheregood]) / ngood)+"\n")

    # accumulate onto net chi2
    allq_X2 += X2
    allq_delta += delta
    allq_badGridCell = np.logical_or(allq_badGridCell, badGridCell)
    ftxt.write("\n")
    ftxt.flush()

ftxt.write("============================\n")
ftxt.write("= Total over all variables =\n")
ftxt.write("============================\n")
# clear unneeded variables
del(data,theory,data2,theory2,stddev2,this_theory,this_data)

# report quality
ftxt.write("shape: "+str(np.shape(allq_X2))+"\n")
allq_wherebad  = np.where(allq_badGridCell)
allq_wheregood = np.where(np.logical_not(allq_badGridCell))
allq_ngood = np.shape(allq_wheregood)[1]
ftxt.write(str(len(allq_wherebad[0]))+" bad points\n")

# sum over cell to get scalar
ftxt.write("total X2/dof: "   +str(np.sum(   allq_X2[allq_wheregood]) / (allq_ngood*len(qlist)))+"\n")
ftxt.write("total delta/dof: "+str(np.sum(allq_delta[allq_wheregood]) / (allq_ngood*len(qlist)))+"\n")
ftxt.write("min fluxfac = "+str(np.min(fluxfac[allq_wheregood]))+"\n")
ftxt.write("max fluxfac = "+str(np.max(fluxfac[allq_wheregood]))+"\n")


##################
# OUTPUT TO HDF5 #
##################
# output to hdf5
ftxt.write("\n")
ftxt.write("writing to "+outfilename+".h5\n")
f = h5py.File(outfilename+".h5", "w")
 
# grids
f.create_dataset("xgrid",data=xgrid)
f.create_dataset("ygrid",data=ygrid)
f.create_dataset("zgrid",data=zgrid)
f.create_dataset("Egrid",data=Egrid)

# for color plots of chi squared
X2_sgAvg = np.sum(allq_X2,axis=(0,4))/(ns*ng*len(qlist))
f.create_dataset("chi2_yz",data=X2_sgAvg[0,:,:])
f.create_dataset("chi2_xz",data=X2_sgAvg[:,0,:])
f.create_dataset("chi2_xy",data=X2_sgAvg[:,:,0])

# for color plots of delta
delta_sgAvg = np.sum(allq_delta,axis=(0,4))/(ns*ng*len(qlist))
f.create_dataset("err_yz",data=delta_sgAvg[0,:,:])
f.create_dataset("err_xz",data=delta_sgAvg[:,0,:])
f.create_dataset("err_xy",data=delta_sgAvg[:,:,0])

# for line plots of delta with energy and species
X2_sg    = np.zeros((ns,ng))
delta_sg = np.zeros((ns,ng))
ngood_sg = np.zeros((ns,ng))
for s in range(ns):
    for ig in range(ng):
        with np.errstate(invalid='ignore'):
            wheregood = np.where(allq_delta[:,:,:,:,ig]>0)
        ngood_sg[s,ig] = np.shape(wheregood)[1]*len(qlist)
        delta_sg[s,ig] = np.sum(allq_delta[:,:,:,:,ig][wheregood])/ngood_sg[s,ig]
        X2_sg[s,ig]    = np.sum(   allq_X2[:,:,:,:,ig][wheregood])/ngood_sg[s,ig]
f.create_dataset("err_sg",data=delta_sg)
f.create_dataset("X2_sg",data=delta_sg)
f.create_dataset("ngood_sg",data=ngood_sg)

# for color grid plots of delta, chi2 with flux factor
fgrid = np.linspace(0,1,nFluxFacBins+1)
f.create_dataset("fluxfacgrid",data=fgrid)
X2_f = np.zeros(nFluxFacBins)
delta_f = np.zeros(nFluxFacBins)
ngood_f = np.zeros(nFluxFacBins)
for fi in range(nFluxFacBins):
    fmin = fgrid[fi]
    fmax = fgrid[fi+1]
    with np.errstate(invalid='ignore'):
        wheregood = np.where((fluxfac>=fmin) & (fluxfac<=fmax) & np.logical_not(allq_badGridCell))
    ngood_f[fi] = np.shape(wheregood)[1]*len(qlist)
    delta_f[fi] = np.sum(allq_delta[wheregood])/ngood_f[fi]
    X2_f[fi]    = np.sum(   allq_X2[wheregood])/ngood_f[fi]
f.create_dataset("X2_fluxfac", data=X2_f)
f.create_dataset("err_fluxfac", data=delta_f)
f.create_dataset("ngood_fluxfac", data=ngood_f)

f.close()

