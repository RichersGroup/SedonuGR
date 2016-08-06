import matplotlib.pyplot as plt
import numpy as np
import math
import os
import sys
import itertools
from pylab import *
from matplotlib.ticker import AutoMinorLocator
import h5py

sigma = 2000

# read in data
filename = "fluid_00001.h5"
f = h5py.File(filename, "r")
nue_distribution = np.array(f.get("distribution(erg|ccm,lab)")[:,0,0,:,:]) # r, mu, phi
nr = len(nue_distribution)
nmu = len(nue_distribution[0])
nphi = len(nue_distribution[0][0])

# set up the grid
rgrid = f.get("r(cm)")
rcenter = np.zeros(nr)
rcenter2 = np.zeros(nr-1)
rcenter3 = np.zeros(nr-2)
for i in range(nr):
    rcenter[i] = 0.5 * (rgrid[i] + rgrid[i+1])
for i in range(nr-1):
    rcenter2[i] = 0.5 * (rcenter[i] + rcenter[i+1])
for i in range(nr-2):
    rcenter3[i] = 0.5 * (rcenter2[i] + rcenter2[i+1])

mugrid = f.get("distribution_costheta_grid(lab)")
mucenter = np.zeros(nmu)
for i in range(nmu):
    mucenter[i] = 0.5 * (mugrid[i] + mugrid[i+1])

# calculate the energy density and flux arrays
edens = np.zeros(nr)
flux = np.zeros(nr)
for i in range(nr):
    for mu in range(nmu):
        for phi in range(nphi):
            edens[i] += nue_distribution[i][mu][phi]
            flux[i]  += nue_distribution[i][mu][phi] * mucenter[mu]

# calculate the gradient of the energy density
r2ddr = np.zeros(nr-1)
gradE = np.zeros(nr-2)
for i in range(nr-1):
    r = rcenter2[i]
    dr     = rcenter[i+1] - rcenter[i]
    dedens =   edens[i+1] -   edens[i]
    r2ddr[i] = r**2 * dedens / dr
for i in range(nr-2):
    r = rcenter3[i]
    dr     = rcenter2[i+1] - rcenter2[i]
    dderiv =    r2ddr[i+1] -    r2ddr[i]
    gradE[i] = 1./r * dderiv / dr

# evaluate Fick's law
fick_flux = -gradE / (3.0 * sigma)

## Set up plots
panelheight = 0.98  # height of a single panel
panelheighttop = 0.7  # height of a single panel
left, bottom, width, height = 0.15, 0.15, 0.8, 0.8
rect0 = [left, bottom, width, height*panelheight]
fig = plt.figure()#figsize=(11., 12.))
font = 20

#rcParams['xtick.major.size'] = 13
#rcParams['xtick.major.pad']  = 8
#rcParams['xtick.minor.size'] = 6
#rcParams['ytick.major.size'] = 13
#rcParams['ytick.major.pad']  = 8
#rcParams['ytick.minor.size'] = 6

# annotation text
fig.gca().get_xaxis().set_visible(False)
fig.gca().get_yaxis().set_visible(False)
plt.setp(fig.gca(),frame_on=False)

min1 = min(flux)
min2 = min(fick_flux)
max1 = max(flux)
max2 = max(fick_flux)
ax0 = fig.add_axes(rect0)
ax0.axis([min(rgrid), max(rgrid), -1e20, 1e20])#min(min1,min2), max(max1,max2)])

minorLocator = AutoMinorLocator()

# draw plots
plist = []
p, = ax0.plot(rcenter,flux, color="k", linewidth=1)
plist.append(p)
p, = ax0.plot(rcenter3,fick_flux, color="b", linewidth=1)
plist.append(p)

#ax0.xaxis.set_minor_locator(minorLocator)

#rc('legend', fontsize=font)

ax0.legend(plist, ["$\mathrm{Flux}\,[\mathrm{erg/ccm}]$",r"$-\nabla E / 3\sigma\,[\mathrm{erg/ccm}]$"],
           loc=(0.5,0.8),).draw_frame(0)

#labels = getp(ax0, 'xticklabels')
#setp(labels, size=font)
#labels = getp(ax0, 'yticklabels')
#setp(labels, size=font)

## For some reason, this has to go after the plots
xlabelx = 0.5
xlabely = -0.15
ylabelx = -0.2
ylabely = 0.5

xlabelstring = '$r\,[\mathrm{cm}]$'
plt.text(xlabelx,xlabely,xlabelstring,fontsize=font,
         horizontalalignment="center",transform=ax0.transAxes)
ylabelstring = r'Energy Density [erg/ccm]'
plt.text(ylabelx,ylabely,ylabelstring,fontsize=font,
         verticalalignment="center",rotation='vertical',transform=ax0.transAxes)


plt.savefig('ficks_law.pdf',bbox_inches='tight')
