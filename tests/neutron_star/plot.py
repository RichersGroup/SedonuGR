import h5py
import numpy as np
import matplotlib.pyplot as plt
k_MeV = 8.6173303e-11 #MeV/K

plt.xlabel("Radius (km)")
plt.ylabel(r"$Y_e$")
for i in range(20):
    # get Sedonu data
    fileid = i+1
    f = h5py.File("fluid_"+str(fileid).zfill(5)+".h5","r")
    r = np.array(f["axes/x0(cm)[mid]"])/1e5
    ye = np.array(f["Ye"])
    plt.plot(r,ye)    
plt.savefig("r_ye.pdf")

plt.xlabel("Radius (km)")
plt.ylabel(r"$T$ (MeV)")
for i in range(20):
    # get Sedonu data
    fileid = i+1
    f = h5py.File("fluid_"+str(fileid).zfill(5)+".h5","r")
    r = np.array(f["axes/x0(cm)[mid]"])/1e5
    T = np.array(f["T_gas(MeV,tet)"])*k_MeV
    plt.plot(r,T)    
plt.savefig("r_T.pdf")
