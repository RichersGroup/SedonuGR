import numpy as np
import h5py
import matplotlib.pyplot as plt
from matplotlib import gridspec
import os
import sys
sys.path.append("scripts/")
from interpolate_table import *

#os.system("~/software/sedonu/sedonu param.lua")


# physical constants
kB = 1.380658e-16 # erg/K
MeV = 1.60218e-6 # erg
h = 6.6260755e-27 # erg s
c = 2.99792458e10 # cm/s

# Fermi-Dirac function
def FD_theory(hnu, mu, kT): # hnu is energy in MeV, mu is chem. pot. in MeV, T is temperature in MeV
    return 1. / (1. + np.exp( (hnu-mu) / kT ))

def FD_data(nu_mid, dnu3_3, edens):
    return edens * c**3 / (4.*np.pi * h*nu_mid * dnu3_3)

def get_data_3D(infile, ix):
    edens0 = np.array(infile["distribution0(erg|ccm,tet)"])[ix[0],ix[1],ix[2],:,0]
    edens1 = np.array(infile["distribution1(erg|ccm,tet)"])[ix[0],ix[1],ix[2],:,0]
    edens2 = np.array(infile["distribution2(erg|ccm,tet)"])[ix[0],ix[1],ix[2],:,0]/4.
    T = np.array(infile["T_gas(K,tet)"])[ix[0],ix[1],ix[2]]*kB / MeV
    rho = np.array(infile["rho(g|ccm,tet)"])[ix[0],ix[1],ix[2]]
    Ye = np.array(infile["Ye"][ix[0],ix[1],ix[2]])
    x = np.array(infile["axes/x0(cm)[edge]"])
    dx = x[-1]-x[-2]
    return edens0, edens1, edens2, rho, T, Ye, dx
def get_data_0D(infile):
    tstep = 0.01
    edens0 = np.array(infile["distribution0(erg|ccm,tet)"])[:,0]/tstep
    edens1 = np.array(infile["distribution1(erg|ccm,tet)"])[:,0]/tstep
    edens2 = np.array(infile["distribution2(erg|ccm,tet)"])[:,0]/(4.*tstep)
    T = np.array(infile["T_gas(K,tet)"])*kB / MeV
    rho = np.array(infile["rho(g|ccm,tet)"])
    Ye = np.array(infile["Ye"])
    dx = 1
    return edens0, edens1, edens2, rho, T, Ye, dx

infile = h5py.File("fluid_00001.h5","r")
nu_edge = np.array(infile["axes/frequency(Hz)[edge]"])
nu_mid = np.array(infile["axes/frequency(Hz)[mid]"])
if("axes/x0(cm)[mid]" in infile):
    edens0,edens1,edens2,rho,T,Ye,dx = get_data_3D(infile,[1,1,1])
else:
    edens0,edens1,edens2,rho,T,Ye,dx = get_data_0D(infile)

dnu = np.array([ (nu_edge[i+1]-nu_edge[i]) for i in range(len(nu_mid)) ])
dnu3_3 = np.array([ (nu_edge[i+1]**3-nu_edge[i]**3)/3. for i in range(len(nu_mid)) ])
print("rho=",rho)
print("T=", T)
print("Ye=",Ye)
infile.close()

# GET EOS DATA
for line in open("inelastic_scatter_kernel.lua"):
    if "nulib_table" in line:
        nulib_filename = line.split("\"")[1]
    if "randomwalk_min_optical_depth" in line:
        rw_tau = float(line.split("=")[1])
    if "nulib_eos" in line:
        eosfilename = line.split("\"")[1]

eosfile = h5py.File(eosfilename,"r")
mue = interpolate_eos(rho, T, Ye, eosfile, "mu_e")
mun = interpolate_eos(rho, T, Ye, eosfile, "mu_n")
mup = interpolate_eos(rho, T, Ye, eosfile, "mu_p")
muhat = interpolate_eos(rho, T, Ye, eosfile, "muhat")
munue = mue - muhat
#munue = -mup #mue - muhat
print("mue=",mue)
print("mup=",mup)
print("mun=",mun)
print("munue=",munue)


# GET NULIB DATA
nulibfile = h5py.File(nulib_filename,"r")
absopac = np.array([[interpolate_eas(ig, s, rho, T, Ye, nulibfile, "absorption_opacity") for ig in range(len(nu_mid))] for s in range(3)])
scatopac = np.array([[interpolate_eas(ig, s, rho, T, Ye, nulibfile, "scattering_opacity") for ig in range(len(nu_mid))] for s in range(3)])
emis = np.array([[interpolate_eas(ig, s, rho, T, Ye, nulibfile, "emissivities") for ig in range(len(nu_mid))] for s in range(3)])
# NuLib emissivity is erg/ccm/s/sr/MeV
NulibFD = [emis[s]*c**2 * (h*dnu/MeV) / (absopac[s] * h*nu_mid**3 * dnu) for s in range(3)]
NulibFD[2] /= 4.
E = h*nu_mid/MeV

# GENERATE THE PLOTS
fig = plt.figure(figsize=(8,8))
gs = gridspec.GridSpec(2, 1)
axes = [plt.subplot(gsi) for gsi in gs]

axes[0].loglog(E, FD_theory(E, munue, T), color="blue", label=r"BB $\nu_e$")
axes[0].loglog(E, FD_theory(E, -munue, T), color="red", label=r"BB $\bar{\nu}_e$")
axes[0].loglog(E, FD_theory(E, 0, T), color="green", label=r"BB $\nu_x$")
axes[0].scatter(E, FD_data(nu_mid, dnu3_3, edens0/1e-1), label=r"Data", color="blue")
axes[0].scatter(E, FD_data(nu_mid, dnu3_3, edens1/1e-1), color="red")
axes[0].scatter(E, FD_data(nu_mid, dnu3_3, edens2/1e-1), color="green")
axes[0].plot(E, NulibFD[0], label=r"NuLib", color="blue", linestyle="--")
axes[0].plot(E, NulibFD[1], color="red", linestyle="--")
axes[0].plot(E, NulibFD[2], color="green", linestyle="--")
axes[0].legend(loc=3,frameon=False)
axes[0].axhline(1,color="gray")
axes[0].set_ylabel(r"f")
axes[0].set_ylim(10**-14,10**1)

axes[1].loglog(E, absopac[0], label=r"$\kappa_a$", color="blue")
axes[1].loglog(E, absopac[1], color="red")
axes[1].loglog(E, absopac[2], color="green")
axes[1].loglog(E, scatopac[0], linestyle="--", linewidth=2, label=r"$\kappa_s$", color="blue")
axes[1].loglog(E, scatopac[1], linestyle="--", linewidth=2, color="red")
axes[1].loglog(E, scatopac[2], linestyle="--", linewidth=2, color="green")
axes[1].axhline(1./dx,color="gray",label=r"$1/dx$")
axes[1].axhline(rw_tau/dx,color="gray",label=r"$\tau_\mathrm{RW}/dx$")
axes[1].legend(loc=5,frameon=False)

axes[1].set_xlabel(r"$h\nu\,(\mathrm{MeV})$")
axes[0].set_title(r"$\log\rho=$%.1f"%(np.log10(rho))+" $T=$%.1f"%T+" $Ye=$%.2f"%Ye+r" $\mu_{\nu_e}=$%.1f"%munue+r" $\mu_e=$%.1f"%mue+r" $\mu_n=$%.1f"%mun+r" $\mu_p=$%.1f"%mup)

for ax in axes:
    ax.set_xlim(.9,400)
plt.subplots_adjust(wspace=0, hspace=0)
plt.savefig("FD_compare.pdf",bbox_inches="tight")
