import matplotlib as mpl
mpl.use('Agg')
import h5py
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, '../')
import tools
import matplotlib.gridspec as gridspec

tolerance = 6e-2
oldfile = h5py.File("4timesHigh_1D_withPairBrems.h5","r")

dist = np.array(oldfile["distribution(erg|ccm,lab)"])
fold = np.array([
    dist[:,0,:,:],
    dist[:,1,:,:],
    dist[:,2,:,:]])
#Hold = np.array(oldfile["H_abs(erg|s|g,com)"]) - np.array(oldfile["C_emit(erg|s|g,com)"])
#dYedt_old = np.array(oldfile["dYe_dt_abs(1|s,lab)"])

newfile = h5py.File("fluid_00001.h5","r")
r = np.array(newfile["axes/x0(cm)[mid]"])
fnew = np.array([
    newfile["distribution0(erg|ccm,tet)"],
    newfile["distribution1(erg|ccm,tet)"],
    newfile["distribution2(erg|ccm,tet)"]])

def normalized_error(string,old, new):
    npoints = np.prod(np.shape(old))
    error = np.sum(new-old) / npoints
    print(string, npoints, error)

    if(np.abs(error)<tolerance):
        return True
    else:
        return False



fig = plt.figure(figsize=(8,12))
gs = gridspec.GridSpec(3, 2)
axes = [plt.subplot(gsi) for gsi in gs]
passing = True

edens_new = fnew[:,:,:,0].sum(axis=(2))
edens_old = fold[:,:,:,0].sum(axis=(2))
edens_old[2] *= 4.
axes[0].loglog(r, edens_old[0],'b-',label=r"$\nu_e$",linewidth=2)
axes[0].loglog(r, edens_old[1],'r-',label=r"$\bar{\nu}_e$", linewidth=2)
axes[0].loglog(r, edens_old[2],'g-',label=r"$\nu_x$",linewidth=2)
axes[0].loglog(r, edens_new[0],'k-',label="NEW")
axes[0].loglog(r, edens_new[1],'k-')
axes[0].loglog(r, edens_new[2],'k-')
axes[0].set_ylabel("Energy Density (erg/ccm)")
axes[0].legend()
passing = passing and normalized_error("edens:",edens_old/edens_old, edens_new/edens_old)

flux_new = fnew[:,:,:,1].sum(axis=(2))
flux_old = fold[:,:,:,1].sum(axis=(2))
flux_old[2] *= 4.
fluxfac_old = flux_old / edens_old
fluxfac_new = flux_new / edens_new
axes[2].semilogx(r, fluxfac_old[0],'b-',linewidth=2)
axes[2].semilogx(r, fluxfac_old[1],'r-',linewidth=2)
axes[2].semilogx(r, fluxfac_old[2],'g-',linewidth=2)
axes[2].semilogx(r, fluxfac_new[0],'b-')
axes[2].semilogx(r, fluxfac_new[1],'r-')
axes[2].semilogx(r, fluxfac_new[2],'g-')
axes[2].set_ylabel("Flux Factor")
passing = passing and normalized_error("fluxfac:",fluxfac_old, fluxfac_new)

prr_new = fnew[:,:,:,2].sum(axis=(2))
prr_old = fold[:,:,:,2].sum(axis=(2))
prr_old[2] *= 4.
prrfac_old = prr_old / edens_old
prrfac_new = prr_new / edens_new
axes[4].semilogx(r, prrfac_old[0],'b-',linewidth=2)
axes[4].semilogx(r, prrfac_old[1],'r-',linewidth=2)
axes[4].semilogx(r, prrfac_old[2],'g-',linewidth=2)
axes[4].semilogx(r, prrfac_new[0],'b-')
axes[4].semilogx(r, prrfac_new[1],'r-')
axes[4].semilogx(r, prrfac_new[2],'g-')
axes[4].set_ylabel("Eddington Factor")
passing = passing and normalized_error("Prrfac:",prrfac_old, prrfac_new)

nugrid = newfile["axes/frequency(Hz)[mid]"]
ndens_new = (fnew[:,:,:,0] / nugrid).sum(axis=(2))
ndens_old = (fold[:,:,:,0] / nugrid).sum(axis=(2))
ndens_old[2] *= 4.
avge_old = edens_old / ndens_old * tools.h / tools.MeV
avge_new = edens_new / ndens_new * tools.h / tools.MeV
axes[1].semilogx(r, avge_old[0],'b-',linewidth=2)
axes[1].semilogx(r, avge_old[1],'r-',linewidth=2)
axes[1].semilogx(r, avge_old[2],'g-',linewidth=2)
axes[1].semilogx(r, avge_new[0],'k-')
axes[1].semilogx(r, avge_new[1],'k-')
axes[1].semilogx(r, avge_new[2],'k-')
axes[1].set_ylabel("Average Energy (MeV)")
passing = passing and normalized_error("avge:",avge_old/avge_old, avge_new/avge_old)

rho = np.array(newfile["rho(g|ccm,tet)"])
Hnet_new = (np.array(newfile["four-force[abs](erg|ccm|s,tet)"])[:,3] + \
           np.array(newfile["four-force[emit](erg|ccm|s,tet)"])[:,3]) / rho
Hnet_old = (np.array(oldfile["H_abs(erg|s|g,com)"]) - np.array(oldfile["C_emit(erg|s|g,com)"]))
Hnet_new[np.where(r<60e5)] = 0
Hnet_old[np.where(r<60e5)] = 0
axes[3].semilogx(r, Hnet_old, 'r-',linewidth=2)
axes[3].semilogx(r, Hnet_new, 'k-')
axes[3].set_ylabel("Net Heating (erg/g/s)")
maxval = np.max(np.abs(Hnet_old))
passing = passing and normalized_error("Hnet",Hnet_old/maxval, Hnet_new/maxval)

Lnet_new = np.array(newfile["l_abs(1|s|ccm,tet)"])
dYedt_new = Lnet_new / (rho / tools.m_n)
dYedt_old = np.array(oldfile["dYe_dt_abs(1|s,lab)"])
axes[5].semilogx(r, dYedt_old, 'r-',linewidth=2)
axes[5].semilogx(r, dYedt_new, 'k-')
axes[5].set_ylabel(r"$dY_e/dt$ from absorption")
maxval = np.max(np.abs(dYedt_old))
passing = passing and normalized_error("dYedt_abs",dYedt_old/maxval, dYedt_new/maxval)

axes[0].xaxis.set_ticklabels([])
axes[1].xaxis.set_ticklabels([])
axes[2].xaxis.set_ticklabels([])
axes[3].xaxis.set_ticklabels([])
axes[1].yaxis.tick_right()
axes[3].yaxis.tick_right()
axes[5].yaxis.tick_right()
axes[1].yaxis.set_label_position("right")
axes[3].yaxis.set_label_position("right")
axes[5].yaxis.set_label_position("right")
axes[4].set_xlabel("Radius (cm)")
axes[5].set_xlabel("Radius (cm)")

for ax in axes:
    ax.axvline(168*1e5,color="gray")
    #ax.set_xlim(5e5,5e6)
plt.subplots_adjust(wspace=0, hspace=0)
plt.savefig("compare.pdf",bbox_inches="tight")

if passing:
    print("SUCCESS")
else:
    raise Exception("one of the above errors is above the tolerance of "+str(tolerance))
