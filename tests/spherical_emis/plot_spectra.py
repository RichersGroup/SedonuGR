import h5py
import numpy as np
import matplotlib.pyplot as plt

tolerance = 0.05
k_b = 1.3806488e-16 # erg/K
c = 2.99e10         #cm/s
h = 6.62606957e-27  #erg/Hz
pi = 3.14159265359
k_MeV = 1.16046e10
MeV = 0.0000016021773
r = 1.5 #cm
rout = 100 #cm
T = 10*k_MeV #K
mue = 10*MeV
rs = 1.0
alpha_core = np.sqrt(1.-rs/r)
alpha_out = np.sqrt(1. - rs/rout)

f = h5py.File("fluid_00001.h5","r")
nu_grid = np.array(f["axes/frequency(Hz)[mid]"])
nu_edge = np.array(f["axes/frequency(Hz)[edge]"])
nu_delta = [nu_edge[i+1]-nu_edge[i] for i in range(nu_grid.size)]
data0 = np.array(f["spectrum0(erg|s)"][:,0,0])/nu_delta/(4.*np.pi)
data1 = np.array(f["spectrum1(erg|s)"][:,0,0])/nu_delta/(4.*np.pi)

def theory(x,mu):
    return pi*r*r*x*x*x*h/c/c*1/(np.exp((h*x-mu)/(k_b*T))+1.)
xgrid = nu_grid
xgrid_shifted = xgrid * alpha_out/alpha_core
theory0 = theory(xgrid, mue)
theory0GR = theory(xgrid_shifted, mue)
theory1 = theory(xgrid, -mue)
theory1GR = theory(xgrid_shifted, -mue)

plt.title(r"$R_\mathrm{core}=$"+str(r)+r"cm   $T=$"+str(T/k_MeV)+r"MeV   $\mu_{\nu_e}=$"+str(mue/MeV)+"MeV")
plt.xlabel("Neutrino Frequency (Hz) (2.5e20 Hz/MeV)")
plt.ylabel(r"$\nu_e$ Energy Flux (erg/s/Hz/sr)")
plt.plot(xgrid  , theory0,   'g--', label=r"$\nu_e$ Newtonian")
plt.plot(xgrid  , theory0GR, 'g',   label=r"$\nu_e$ GR")
plt.plot(nu_grid, data0,     'go',  label=r"$\nu_e$ Sedonu")
plt.plot(xgrid  , theory1,   'b--', label=r"$\bar{\nu}_e$ Newtonian")
plt.plot(xgrid  , theory1GR, 'b',   label=r"$\bar{\nu}_e$ GR")
plt.plot(nu_grid, data1,     'bo',  label=r"$\bar{\nu}_e$ Sedonu")
plt.legend()
plt.savefig("compare_spectra.pdf")

# do error checking
with open("output.txt") as search:
    for line in search:
        if "DO_GR" in line:
            do_gr = int(line[-2])
if(do_gr==0):
    error0 = np.sum(np.abs(data0-theory0)) / np.sum(data0+theory0)
    error1 = np.sum(np.abs(data1-theory1)) / np.sum(data1+theory1)
else:
    error0 = np.sum(np.abs(data0-theory0GR)) / np.sum(data0+theory0GR)
    error1 = np.sum(np.abs(data1-theory1GR)) / np.sum(data1+theory1GR)
print("error for species 0 =",error0)
print("error for species 1 =",error1)
if error0>tolerance or error1>tolerance:
    raise Exception("spherical_emis results are outside of the tolerance.")
else:
    print("SUCCESS")
