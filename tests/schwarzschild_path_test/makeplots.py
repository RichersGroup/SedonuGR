import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.ticker as ticker

mpl.rcParams['font.size'] = 18
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['xtick.major.size'] = 7
mpl.rcParams['xtick.major.width'] = 2
mpl.rcParams['xtick.major.pad'] = 8
mpl.rcParams['xtick.minor.size'] = 4
mpl.rcParams['xtick.minor.width'] = 2
mpl.rcParams['ytick.major.size'] = 7
mpl.rcParams['ytick.major.width'] = 2
mpl.rcParams['ytick.minor.size'] = 4
mpl.rcParams['ytick.minor.width'] = 2
mpl.rcParams['axes.linewidth'] = 2

fig = plt.figure()

##############
# AROUND.PDF #
##############

fig, axes = plt.subplots(2, sharex=True, gridspec_kw={'hspace':0, 'height_ratios':[2, 1]},figsize=(6,6))
xarr = np.arange(0,1,.05)
axes[0].plot(xarr, 1.5*np.ones(xarr.shape),color='k',linestyle="--")
axes[1].plot(xarr, np.ones(xarr.shape),color='k',linestyle="--")
#axes[0].fill_between(np.arange(0,1,.05),0,1,color='black')
#axes[0].text(.12,.9,"Black Hole", color='white',ha='center', va='center',fontsize=18)
#axes[0].set_ylim(1.484,1.506)
#axes[1].set_ylim(.85,1.04)
axes[0].yaxis.set_major_locator(ticker.MultipleLocator(0.01))
axes[0].set_xlim(0,0.25)
axes[0].set_ylabel(r"$r/r_\mathrm{sch}$")
axes[1].set_ylabel(r"$k_t/k_{t0}$")
axes[1].set_xlabel(r"Orbits")
for ax in axes:
    ax.minorticks_on()
    ax.tick_params(axis='both',which='both',direction='in',top=True,right=True)

tag = ['1','2','4','_smallstep2','_smallstep4'] #[72,36,18,9]
colors=['k','g','r','b','orange','pink']
linestyles = ['-','-','-','-','-','-']
markers = [".",'.','.','.','.']
plist = []
key = ["Base",r"Res. $\times 2$", r"Res. $\times 4$",r'Res. $\times 4$ Step Size $\times 0.5$',r'Res. $\times 4$ Step Size $\times 0.25$']
for i in range(len(tag)):
    filename = "around"+tag[i]+".dat"
    x = np.genfromtxt(filename, usecols=(0))
    y = np.genfromtxt(filename, usecols=(1))
    r = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y,x)
    orbits = phi/(2.*np.pi)
    p, = axes[0].plot(orbits, r, color=colors[i], linestyle=linestyles[i],marker=markers[i])
    plist.append(p)

    klowt = np.genfromtxt(filename, usecols=(17))
    axes[1].plot(orbits, klowt/klowt[0], color=colors[i], linestyle=linestyles[i],marker=markers[i])

axes[0].legend(plist,key,frameon=False,loc=3,ncol=1,fontsize=16)


plt.savefig("around.pdf", bbox_inches='tight')

##############
# RADIAL.PDF #
##############
plt.clf()

fig, axes = plt.subplots(2, sharex=True, gridspec_kw={'hspace':0, 'height_ratios':[2, 1]},figsize=(6,6))
rarr = np.arange(1.5,10,.02)
sqrtr = np.sqrt(rarr)
tarr = rarr + np.log(np.abs(rarr-1.))# from wikipedia page on geodesics
tarr = tarr - tarr[0]


axes[0].plot(tarr, np.zeros(tarr.shape),color='k',linestyle="--")
axes[1].plot(tarr, np.ones(tarr.shape),color='k',linestyle="--")
#axes[0].fill_between(np.arange(0,1,.05),0,1,color='black')
#axes[0].text(.12,.9,"Black Hole", color='white',ha='center', va='center',fontsize=18)
axes[0].set_ylim(-.25,0.02)
axes[1].set_ylim(.98,1.039)
axes[0].set_xlim(0,10)
axes[0].set_ylabel(r"$(r-r_\mathrm{analytic})/r_\mathrm{sch}$")
axes[1].set_ylabel(r"$k_t/k_{t0}$")
axes[1].set_xlabel(r"$ct/r_\mathrm{sch}$")
for ax in axes:
    ax.minorticks_on()
    ax.tick_params(axis='both',which='both',direction='in',top=True,right=True)

tag = ['1','2','4','_smallstep2','_smallstep4'] #[72,36,18,9]
colors=['k','g','r','b','orange','pink']
linestyles = ['-','-','-','-','-','-']
markers = ['.','.','.','.','']
plist = []
key = ["Base",r"Res. $\times 2$", r"Res. $\times 4$",r'Res. $\times 4$ Step Size $\times 0.5$',r'Res. $\times 4$ Step Size $\times 0.25$']
for i in range(len(tag)):
    filename = "radial"+tag[i]+".dat"
    x = np.genfromtxt(filename, usecols=(0))
    y = np.genfromtxt(filename, usecols=(1))
    t = np.genfromtxt(filename, usecols=(3))
    r = np.sqrt(x**2 + y**2)
    p, = axes[0].plot(t, r-np.interp(t, tarr,rarr), color=colors[i], linestyle=linestyles[i],marker=markers[i])
    plist.append(p)

    klowt = np.genfromtxt(filename, usecols=(17))
    axes[1].plot(t, klowt/klowt[0], color=colors[i], linestyle=linestyles[i],marker=markers[i])

rgrid = np.genfromtxt("empty_sphere1.mod", usecols=(0))
for r in rgrid:
    for ax in axes:
        ax.axvline(np.interp(r, rarr,tarr),linewidth=0.5, color="gray")

axes[0].text(6.3,-.1,r"$t_\mathrm{analytic}$ when",fontsize=10,rotation=-90,horizontalalignment='center',verticalalignment='center',color="gray")
axes[0].text(5.9,-.1,r"$r_\mathrm{analytic} = r_\mathrm{out}$",fontsize=10,rotation=-90,horizontalalignment='center',verticalalignment='center',color="gray")

#plist[3],plist[4] = plist[4],plist[3]
#key[3],key[4] = key[4],key[3]
axes[0].legend(plist,key,frameon=False,loc=(.04,0.04),ncol=1,fontsize=14)

plt.savefig("radial.pdf", bbox_inches='tight')
