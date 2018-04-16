import h5py
import numpy as np
import matplotlib.pyplot as plt

# get Sedonu data
f = h5py.File("fluid_00001.h5","r")
r = np.array(f["axes/x0(cm)[mid]"])/1e5
ye = np.array(f["Ye"])
rho = np.array(f["rho(g|ccm,tet)"])

# get EOS data
eos_data_filename = "eos.dat"
eos_r = np.genfromtxt(eos_data_filename,usecols=(0))/1e5
eos_rho = np.genfromtxt(eos_data_filename,usecols=(1))
eos_T = np.genfromtxt(eos_data_filename,usecols=(2))
eos_Ye = np.genfromtxt(eos_data_filename,usecols=(3))

plt.xlabel("Radius (km)")
plt.ylabel(r"$Y_e$")
plt.plot(eos_r,eos_Ye)
plt.plot(r,ye,'ro')
plt.savefig("r_ye.pdf")

plt.cla()
plt.xlabel(r"$\rho$ (g/cm)")
plt.ylabel("$Y_e$")
plt.plot(eos_rho,eos_Ye)
plt.plot(rho,ye,'ro')
plt.savefig("rho_ye.pdf")
