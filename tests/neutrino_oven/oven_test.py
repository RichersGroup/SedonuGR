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
axes[0].legend(loc=4)
axes[0].set_ylim(0,1.1*max(np.max(edens),np.max(edens_theory)))
error = edens-edens_theory
error2 = (edens-edens_theory)**2
std_dev = np.sqrt(np.sum(error2)/len(error2) - (np.sum(error)/len(error))**2)
abserror = np.sum(np.abs(error))/len(error)
print("Energy density: error="+str(abserror)+" stddev="+str(std_dev))
if abserror>2.*std_dev:
    isError = True

e_abs = np.array(f["four-force[abs](erg|ccm|s,tet)"])[:,3]
e_emit = np.array(f["four-force[emit](erg|ccm|s,tet)"])[:,3]
dEdt = e_abs + e_emit
axes[1].set_ylabel(r"$dE_\mathrm{int}/dt$ (erg/ccm/s)")
axes[1].plot(r, -e_emit, label="emit")
axes[1].scatter(r, e_abs, label="abs")
axes[1].legend(loc=4)
axes[1].set_ylim(0,1.1*max(np.max(e_abs),np.max(-e_emit)))
error = dEdt
error2 = dEdt**2
std_dev = np.sqrt(np.sum(error2)/len(error2) - (np.sum(error)/len(error))**2)
abserror = np.sum(np.abs(error))/len(error)
print("dEdt: error="+str(abserror)+" stddev="+str(std_dev))
if abserror>2.*std_dev:
    isError = True

l_abs = np.array(f["l_abs(1|s|ccm,tet)"])
l_emit = np.array(f["l_emit(1|s|ccm,tet)"])
dNdt = l_abs - l_emit
axes[2].set_ylabel(r"$dn_l/dt$ (1/ccm/s)")
axes[2].plot(r, l_emit, label="emit")
axes[2].scatter(r, l_abs, label="abs")
axes[2].legend()
error = dNdt
error2 = dNdt**2
std_dev = np.sqrt(np.sum(error2)/len(error2) - (np.sum(error)/len(error))**2)
abserror = np.sum(np.abs(error))/len(error)
print("dNdt: error="+str(abserror)+" stddev="+str(std_dev))
if abserror>2.*std_dev:
    isError = True

for ax in axes:
    ax.axvline(1.5)
    
axes[2].set_xlabel("$r$ (km)")
axes[0].xaxis.set_ticklabels([])
axes[1].xaxis.set_ticklabels([])

plt.subplots_adjust(wspace=0, hspace=0)
plt.savefig("edens.pdf", bbox_inches="tight")

if isError:
    raise Exception("temperature_redshift results are outside of the tolerance.")
else:
    print("SUCCESS")
