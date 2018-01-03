import h5py
import numpy as np
import matplotlib.pyplot as plt

k_b = 1.3806488e-16 # erg/K
c = 2.99e10         #cm/s
h = 6.62606957e-27  #erg/Hz
pi = 3.14159265359
k_MeV = 1.16046e10
MeV = 0.0000016021773
r = 9.99999e5                             #cm
T = 10*k_MeV                        #K
mue = 10*MeV

f = h5py.File("fluid_00001.h5","r")
nu_grid = np.array(f["nu_grid(Hz)[mid]"])
nu_edge = np.array(f["nu_grid(Hz)[edge]"])
nu_delta = [nu_edge[i+1]-nu_edge[i] for i in range(nu_grid.size)]
data0 = np.array(f["spectrum0(erg|ccm,lab)"][:,0,0])/nu_delta/(4.*np.pi)
data1 = np.array(f["spectrum1(erg|ccm,lab)"][:,0,0])/nu_delta/(4.*np.pi)
data2 = np.array(f["spectrum2(erg|ccm,lab)"][:,0,0])/nu_delta/(4.*np.pi)

def theory(x,mu):
    return pi*r*r*x*x*x*h/c/c*1/(np.exp((h*x-mu)/(k_b*T))+1.)   
xgrid = np.linspace(nu_edge[0], nu_edge[-1], num=200)
theory0 = theory(xgrid, mue)
theory1 = theory(xgrid, -mue)
theory2 = 4.*theory(xgrid, 0)


plt.xlabel("Neutrino Frequency (Hz) (2.5e20 Hz/MeV)")
plt.ylabel("Energy Flux (erg/s/Hz/sr)")
plt.plot(xgrid,theory0)
plt.plot(nu_grid,data0, 'go')
plt.savefig("compare_spectrum_0.pdf")

plt.cla()
plt.xlabel("Neutrino Frequency (Hz) (2.5e20 Hz/MeV)")
plt.ylabel("Energy Flux (erg/s/Hz/sr)")
plt.plot(xgrid,theory1)
plt.plot(nu_grid,data1,'go')
plt.savefig("compare_spectrum_1.pdf")

plt.cla()
plt.xlabel("Neutrino Frequency (Hz) (2.5e20 Hz/MeV)")
plt.ylabel("Energy Flux (erg/s/Hz/sr)")
plt.plot(xgrid,theory2)
plt.plot(nu_grid,data2,'go')
plt.savefig("compare_spectrum_2.pdf")
