#!/usr/bin/python

#import sys
#import subprocess

from numpy import *
from pylab import *
from plot_defaults import *

## get parameters
f = open("eas_T.dat","r")
firstline = f.readline()
pars = firstline.strip('# \r\n')
f.close()

## Open Files
data = loadtxt("eas_T.dat")
outfile = "eas_T.pdf"

## Figure Positioning
fig = figure(figsize=(10,10))
left, bottom, width, height = 0.15, 0.12, 0.82, 0.31
panels=3
heights = [0.9*height, 0.9*height, 0.9*height]
rect, ax = setup_plot_with_vertical_panels(fig, left, bottom, heights, width, panels)


#############
# plot data #
#############
pEmis,     = ax[0].plot(data[:,0],data[:,1],'k',linewidth=3)
pAbsOpac,  = ax[1].plot(data[:,0],data[:,2],'k',linewidth=3)
pScatOpac, = ax[2].plot(data[:,0],data[:,3],'k',linewidth=3)

# Set up axes
xmin = min(data[:,0])
xmax = max(data[:,0])
ax[0].set_xscale('log')
ax[1].set_xscale('log')
ax[2].set_xscale('log')
ax[0].grid(True)
ax[1].grid(True)
ax[2].grid(True)

y0min = min(data[:,1])
y0max = max(data[:,1])
ax[0].set_yscale('log')
ax[0].axis([xmin,xmax,y0min,y0max])
for label in ax[0].get_yticklabels():
    label.set_fontsize(15)

y1min = min(data[:,2])
y1max = max(data[:,2])
ax[1].set_yscale('log')
ax[1].axis([xmin,xmax,y1min,y1max])
for label in ax[1].get_yticklabels():
    label.set_fontsize(15)

y2min = min(data[:,3])
y2max = max(data[:,3])
ax[2].set_yscale('log')
ax[2].axis([xmin,xmax,y2min,y2max])
for label in ax[2].get_yticklabels():
    label.set_fontsize(15)

xmajorFormatter = FormatStrFormatter('%g')
ax[0].xaxis.set_major_formatter(xmajorFormatter)
ymajorFormatter = FormatStrFormatter('%g')
ax[0].yaxis.set_major_formatter(ymajorFormatter)
matplotlib.rc('ytick', labelsize=10)

# Axis Labels
labelsize=30

xlabelx = 0.5
xlabely = -0.3
xlabelstring = r'T (MeV)'
text(xlabelx,xlabely,xlabelstring,fontsize=labelsize,
          horizontalalignment="center",transform=ax[0].transAxes)

y0labelx = -0.15
y0labely = 0.5
y0labelstring = r'$\epsilon$ (erg/cm$^3$/s/ster)'
text(y0labelx,y0labely,y0labelstring,fontsize=0.8*labelsize,
         verticalalignment="center",rotation='vertical',transform=ax[0].transAxes)

y1labelx = -0.15
y1labely = 1.5
y1labelstring = r'$\kappa_a$ (cm$^2$/g)'
text(y1labelx,y1labely,y1labelstring,fontsize=0.8*labelsize,
         verticalalignment="center",rotation='vertical',transform=ax[0].transAxes)

y2labelx = -0.15
y2labely = 2.5
y2labelstring = r'$\kappa_s$ (cm$^2$/g)'
text(y2labelx,y2labely,y2labelstring,fontsize=0.8*labelsize,
         verticalalignment="center",rotation='vertical',transform=ax[0].transAxes)

titlex = 0.5
titley = 3.05
text(titlex,titley,pars,fontsize=0.5*labelsize,
     horizontalalignment="center",transform=ax[0].transAxes)

## Output
plt.savefig(outfile)

