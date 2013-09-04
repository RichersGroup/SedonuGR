#!/usr/bin/python

#import sys
#import subprocess

from numpy import *
from pylab import *
#from matplotlib import *
#from scipy import interpolate
from plot_defaults import *


## Open Files
GRM1r    = loadtxt("../../M1/GRM1/post_processing/profile_rho_0.075.dat")
GRM1y    = loadtxt("../../M1/GRM1/post_processing/profile_ye_0.075.dat")
GRM1s    = loadtxt("../../M1/GRM1/post_processing/profile_temperature_0.075.dat")

NM1r     = loadtxt("../../M1/NM1/post_processing/profile_rho_0.075.dat")
NM1y     = loadtxt("../../M1/NM1/post_processing/profile_ye_0.075.dat")
NM1s     = loadtxt("../../M1/NM1/post_processing/profile_temperature_0.075.dat")

GRLeakr  = loadtxt("../../Leakage/GRLeakage/post_processing/profile_rho_0.075.dat")
GRLeaky  = loadtxt("../../Leakage/GRLeakage/post_processing/profile_ye_0.075.dat")
GRLeaks  = loadtxt("../../Leakage/GRLeakage/post_processing/profile_temperature_0.075.dat")

NLeakr   = loadtxt("../../Leakage/NLeakage/post_processing/profile_rho_0.075.dat")
NLeaky   = loadtxt("../../Leakage/NLeakage/post_processing/profile_ye_0.075.dat")
NLeaks   = loadtxt("../../Leakage/NLeakage/post_processing/profile_temperature_0.075.dat")

GRMBnr   = loadtxt("../../MB/GRMB_no_entropy/post_processing/profile_rho_0.075.dat")
GRMBny   = loadtxt("../../MB/GRMB_no_entropy/post_processing/profile_ye_0.075.dat")
GRMBns   = loadtxt("../../MB/GRMB_no_entropy/post_processing/profile_temperature_0.075.dat")

NMBnr    = loadtxt("../../MB/NMB_no_entropy/post_processing/profile_rho_0.075.dat")
NMBny    = loadtxt("../../MB/NMB_no_entropy/post_processing/profile_ye_0.075.dat")
NMBns    = loadtxt("../../MB/NMB_no_entropy/post_processing/profile_temperature_0.075.dat")

outfile = "../figures/profile.pdf"


## Figure Positioning
fig = figure(figsize=(10,10))
left, bottom, width, height = 0.15, 0.12, 0.82, 0.31
panels=3
heights = [0.9*height, 0.9*height, 0.9*height]
rect, ax = setup_plot_with_vertical_panels(fig, left, bottom, heights, width, panels)


###########
# Density #
###########
pGRM1r,   = ax[0].plot(GRM1r[:,0]/1e5,log10(GRM1r[:,1]),'k',linewidth=3)
pGRLeakr, = ax[0].plot(GRLeakr[:,0]/1e5,log10(GRLeakr[:,1]),'b',linewidth=3)
pGRMBnr,  = ax[0].plot(GRMBnr[:,0]/1e5,log10(GRMBnr[:,1]),'r',linewidth=3)

pNM1r,   = ax[0].plot(NM1r[:,0]/1e5,log10(NM1r[:,1]),'k--',linewidth=3)
pNLeakr, = ax[0].plot(NLeakr[:,0]/1e5,log10(NLeakr[:,1]),'b--',linewidth=3)
pNMBnr,  = ax[0].plot(NMBnr[:,0]/1e5,log10(NMBnr[:,1]),'r--',linewidth=3)


#####################
# Electron Fraction #
#####################
pGRM1y,   = ax[1].plot(GRM1y[:,0]/1e5,GRM1y[:,1],'k',linewidth=3)
pGRLeaky, = ax[1].plot(GRLeaky[:,0]/1e5,GRLeaky[:,1],'b',linewidth=3)
pGRMBny,  = ax[1].plot(GRMBny[:,0]/1e5,GRMBny[:,1],'r',linewidth=3)

pNM1y,   = ax[1].plot(NM1y[:,0]/1e5,NM1y[:,1],'k--',linewidth=3)
pNLeaky, = ax[1].plot(NLeaky[:,0]/1e5,NLeaky[:,1],'b--',linewidth=3)
pNMBny,  = ax[1].plot(NMBny[:,0]/1e5,NMBny[:,1],'r--',linewidth=3)


###############
# Temperature #
###############
pGRM1s,   = ax[2].plot(GRM1s[:,0]/1e5,GRM1s[:,1],'k',linewidth=3)
pGRLeaks, = ax[2].plot(GRLeaks[:,0]/1e5,GRLeaks[:,1],'b',linewidth=3)
pGRMBns,  = ax[2].plot(GRMBns[:,0]/1e5,GRMBns[:,1],'r',linewidth=3)

pNM1s,   = ax[2].plot(NM1s[:,0]/1e5,NM1s[:,1],'k--',linewidth=3)
pNLeaks, = ax[2].plot(NLeaks[:,0]/1e5,NLeaks[:,1],'b--',linewidth=3)
pNMBns,  = ax[2].plot(NMBns[:,0]/1e5,NMBns[:,1],'r--',linewidth=3)


# Set up axes
xmin = 5
xmax = 500
ax[0].set_xscale('log')
ax[1].set_xscale('log')
ax[2].set_xscale('log')

y0min = 7
y0max = 14.99
#ax[0].set_yscale('log')
ax[0].axis([xmin,xmax,y0min,y0max])

y1min = 0
y1max = 0.599
#ax[1].set_yscale('log')
ax[1].axis([xmin,xmax,y1min,y1max])

y2min = 0
y2max = 30
ax[2].axis([xmin,xmax,y2min,y2max])

ax[0].xaxis.set_major_formatter(FuncFormatter(supermongolikeN4a))
ymajorFormatter = FormatStrFormatter('%g')
ax[0].yaxis.set_major_formatter(ymajorFormatter)


# Axis Labels
labelsize=30

xlabelx = 0.5
xlabely = -0.3
xlabelstring = "Radius [km]"
text(xlabelx,xlabely,xlabelstring,fontsize=labelsize,
          horizontalalignment="center",transform=ax[0].transAxes)

y0labelx = -0.17
y0labely = 0.5
y0labelstring = r'$\log_{10}\left(\frac{\rho} {\small\mathrm{g\,cm^{-3}}}\right)$'
text(y0labelx,y0labely,y0labelstring,fontsize=labelsize,
         verticalalignment="center",rotation='vertical',transform=ax[0].transAxes)

y1labelx = -0.15
y1labely = 1.5
y1labelstring = r'$Y_e$'
text(y1labelx,y1labely,y1labelstring,fontsize=labelsize,
         verticalalignment="center",rotation='vertical',transform=ax[0].transAxes)

y2labelx = -0.15
y2labely = 2.5
y2labelstring = r'T $[$MeV$]$'
text(y2labelx,y2labely,y2labelstring,fontsize=labelsize,
         verticalalignment="center",rotation='vertical',transform=ax[0].transAxes)


## Legend
mylegend = ax[2].legend((pGRM1r,pGRLeakr,pGRMBnr,pGRM1r,pNM1r),
                     ('Two-Moment','Leakage','Lightbulb','GR','Newtonian'),
                     loc=(0.32,0.52),axespad=0.05,labelsep=0.01,markerscale=0,ncol=2)
mylegend.draw_frame(False)
mylegendtext = mylegend.get_texts()
mylegendlines = mylegend.get_lines()
setp(mylegendtext, fontsize=20)
#setp(mylegendlines, linewidth=4)
gca().add_artist(mylegend)

## Output
plt.savefig(outfile)

