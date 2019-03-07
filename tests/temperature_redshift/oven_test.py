import matplotlib as mpl
mpl.use('Agg')
import h5py
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, '../')
import tools
import matplotlib.gridspec as gridspec

tolerance = 0.05

with open("param.lua") as search:
    for line in search:
        if "r_core" in line:
            r_core = float(line.split()[2])/1e5
        if("T_core" in line):
            T_core = float(line.split("=")[1].split(",")[0].strip(" ").strip("{"))
        if("core_chem_pot" in line):
            munue = float(line.split("=")[1].split(",")[0].strip(" ").strip("{"))

f = h5py.File("fluid_00001.h5","r")
r = np.array(f["axes/x0(cm)[mid]"])/1e5
T_gas = np.array(f["T_gas(K,tet)"])*tools.k_b/tools.MeV
edens = np.array(f["distribution0(erg|ccm,tet)"]).sum(axis=(1,2,3))
edens_theory = np.array([tools.edens(T_gas[i]*tools.MeV, munue*tools.MeV) for i in range(len(T_gas))])
#edens_theory = np.array([tools.energy_blackbody(T_gas[i]*tools.MeV, munue*tools.MeV) for i in range(len(T_gas))])

fig = plt.figure(figsize=(8,12))
gs = gridspec.GridSpec(3, 1)
axes = [plt.subplot(gsi) for gsi in gs]
axes[0].set_title(r"$R_\mathrm{core}=$"+str(r_core)+r"km   $T_\mathrm{core}=$"+str(T_core)+r"MeV   $\mu_{\nu_e}=$"+str(munue/tools.MeV)+"MeV")

isError = False

axes[0].set_ylabel(r"Energy Density (erg/ccm)")
axes[0].plot(r, edens_theory, "g-",label="theory")
axes[0].scatter(r, edens, label="Sedonu")
error = np.sum(edens-edens_theory)/len(edens)
error2 = np.sum((edens-edens_theory)**2)/len(edens)
std_dev = np.sqrt(error2 - error**2)
print("Energy density: error="+str(error)+" stddev="+str(std_dev))
if np.abs(error)>std_dev:
    isError = True

dEdt = np.array(f["four-force[abs](erg|ccm|s,tet)"])[:,3] + np.array(f["four-force[emit](erg|ccm|s,tet)"])[:,3]
axes[1].set_ylabel(r"$dE_\mathrm{int}/dt$ (erg/ccm/s)")
axes[1].axhline(0,color="g")
axes[1].scatter(r, dEdt)
error = np.sum(dEdt)/len(dEdt)
error2 = np.sum(dEdt**2)/len(dEdt)
std_dev = np.sqrt(error2 - error**2)
print("dEdt: error="+str(error)+" stddev="+str(std_dev))
if np.abs(error)>std_dev:
    isError = True

dNdt = np.array(f["l_abs(1|s|ccm,tet)"]) - np.array(f["l_emit(1|s|ccm,tet)"])
axes[2].set_ylabel(r"$dn_l/dt$ (1/ccm/s)")
axes[2].axhline(0,color="g")
axes[2].scatter(r, dNdt)
error = np.sum(dNdt)/len(dNdt)
error2 = np.sum(dNdt**2)/len(dNdt)
std_dev = np.sqrt(error2 - error**2)
print("dNdt: error="+str(error)+" stddev="+str(std_dev))
if np.abs(error)>std_dev:
    isError = True

axes[2].set_xlabel("$r$ (km)")
axes[0].legend(ncol=2,frameon=False)
axes[0].xaxis.set_ticklabels([])
axes[1].xaxis.set_ticklabels([])

plt.subplots_adjust(wspace=0, hspace=0)
plt.savefig("edens.pdf", bbox_inches="tight")

if isError:
    raise Exception("temperature_redshift results are outside of the tolerance.")
else:
    print("SUCCESS")
