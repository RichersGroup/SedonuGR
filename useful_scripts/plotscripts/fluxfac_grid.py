import matplotlib.pyplot as plt
import h5py
import numpy as np
import os

qsetlist = [
    ['Pxx','Pyy','Pzz'],
    ['Pxy','Pxz','Pyz'],
    ['Lxxx','Lyyy','Lzzz'],
    ['Lxxy','Lxxz','Lxyy','Lxzz','Lxyz','Lyyz','Lyzz']
]

outdir = "plots/f_chi_grid/"
if not os.path.exists(outdir):
    os.makedirs(outdir)

nchi = 20
delta = (2./3.)/nchi
chimid = np.linspace(1./3., 1.0, num=nchi)
chigrid = chimid - delta/2.
chigrid = np.append(chigrid, chimid[-1]+delta/2.)

chistr = []
for chi in chimid:
    chistr.append("%.3f"%chi)

def plot_data(data_good, data_bad,label):
    plt.xlabel("Flux Factor")
    plt.ylabel("Approximation Parameter")
    vmax = max(np.max(data_good),np.max(data_bad))
    plt.pcolormesh(fgrid,chigrid,data_bad,vmin=0,vmax=vmax)#,alpha=0.25)
    plt.pcolormesh(fgrid,chigrid,data_good,vmin=0,vmax=vmax)
    plt.axes().set_aspect('equal', 'datalim')
    plt.colorbar(label=label)
    plt.axes().axis('tight')
    print(np.min(data_good),np.max(data_good))

    fluxfac = np.linspace(0.0,1.0, num=100)

    chi = (1. + 2.*fluxfac**2) / 3.
    plt.plot(fluxfac,chi,'b-',label="Kershaw",linewidth=3)

    chi = 1./3. - 1./3.*fluxfac + fluxfac**2
    plt.plot(fluxfac,chi,'k:',label="Wilson",linewidth=3)

    chi = (3 + 4.*fluxfac**2) / (5. + 2.*np.sqrt(4. - 3.*fluxfac**2))
    plt.plot(fluxfac,chi,'k.-',label="Levermore",linewidth=3)

    chi = 1./3. * (1.0 -2.*fluxfac + 4.*fluxfac**2)
    plt.plot(fluxfac,chi,'k--',label="MEFD",linewidth=3)

    chi = 1./3. + (2.*fluxfac**2)/15. * (3. - fluxfac + 3.*fluxfac**2)
    plt.plot(fluxfac,chi,'k-',label="ME",linewidth=3)

    plt.legend(bbox_to_anchor=(1.35, 1),loc=2)

    
for qset in qsetlist:
    qsetstring = ""
    for q in qset:
        qsetstring = qsetstring+q
    print(qsetstring)
    
    plt.rcParams.update({'font.size': 20})
    
    X2data = []
    errdata = []
    for chi in chistr:
        infilename = "analysis/"+chi+"_s012_g0-16_nf9_"+qsetstring+".h5"
        f = h5py.File(infilename,"r")
        fgrid = f["fluxfacgrid"]
        X2data.append(f["X2_fluxfac"])
        errdata.append(f["err_fluxfac"])
    X2data = np.array(X2data)
    errdata = np.array(errdata)
    errdata_good = np.ma.masked_array(errdata, X2data < 1)
    errdata_bad  = np.ma.masked_array(errdata, X2data >= 1)
    X2data_good  = np.ma.masked_array( X2data, X2data < 1)
    X2data_bad   = np.ma.masked_array( X2data, X2data >= 1)
    # HACK
    X2data[np.where(X2data!=X2data)] = 0
    errdata[np.where(errdata!=errdata)] = 0
    
    plt.clf()
    plot_data(X2data_good,X2data_bad,r"$\chi^2/\nu$")
    plt.savefig(outdir+"/X2_"+qsetstring+".png",bbox_inches="tight")

    plt.clf()
    plot_data(errdata_good,errdata_bad,r"Relative Error")
    plt.savefig(outdir+"/err_"+qsetstring+".png",bbox_inches="tight")
