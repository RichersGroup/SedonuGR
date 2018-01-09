import matplotlib.pyplot as plt
import numpy as np
import math
import os
import sys
from pylab import *

fontsize = 28
fontname = 'times new roman'
font = { 'family' : fontname,
         'weight' : 'normal',
         'size'   : fontsize}
lw=4
bw=3
tw=2
linecolor='r'
pointcolor='k+'
markersize=10
mew=1.5

tol = 1e-3
center_rho = 1e10
center_temp = 3.1623
center_ye = 0.3
epsmin = 1e26
epsmax = 1e35

zR = np.loadtxt('results.dat', skiprows=1, usecols=(0,1,2,4,5,6))
zR = np.transpose(zR)
rhoR   = zR[0] # g/ccm
tempR  = zR[1] # MeV
yeR    = zR[2]
epsR   = zR[3] + zR[4] + zR[5] # erg/ccm

zP = np.loadtxt('predicted.dat', skiprows=1, usecols=(0,1,2,3))
zP = np.transpose(zP)
epsP   = zP[0] # erg/ccm
rhoP   = zP[1] # g/ccm
tempP  = zP[2] # MeV
yeP    = zP[3]

## Set up plots
panelheight = 0.98  # height of a single panel
panelheighttop = 0.7  # height of a single panel
left, bottom, width, height = 0.15, 0.15, 0.27, 0.8
rect0 = [left          , bottom, width, height*panelheight]
rect1 = [left +   width, bottom, width, height*panelheight]
rect2 = [left + 2*width, bottom, width, height*panelheight]
fig = plt.figure(figsize=(22., 5.))

rcParams['xtick.major.size'] = 13
rcParams['xtick.major.pad']  = 8
rcParams['xtick.minor.size'] = 7
rcParams['ytick.major.size'] = 13
rcParams['ytick.major.pad']  = 8
rcParams['ytick.minor.size'] = 7
rcParams['axes.linewidth'] = 3
rcParams['font.family'] = "serif"
rcParams['font.size'] = 32
rcParams['text.usetex'] = True
rcParams['font.serif'] = "palatino"

# annotation text
fig.gca().get_xaxis().set_visible(False)
fig.gca().get_yaxis().set_visible(False)
plt.setp(fig.gca(),frame_on=False)

ax0 = fig.add_axes(rect0)
ax0.set_xscale('log')
ax0.set_yscale('log')
ax0.axis([1e6,1e15,epsmin,epsmax])

# draw plots
indicesR = (np.absolute(tempR-center_temp) < tol) & (np.absolute(yeR-center_ye) < tol)
indicesP = (np.absolute(tempP-center_temp) < tol) & (np.absolute(yeP-center_ye) < tol)
p_Prho, = ax0.plot(rhoP[indicesP],epsP[indicesP],linecolor,linewidth=lw)
p_Rrho, = ax0.plot(rhoR[indicesR],epsR[indicesR],pointcolor,markersize=markersize,mew=mew)

ax0.spines['bottom'].set_linewidth(bw)
ax0.spines['top'].set_linewidth(bw)
ax0.spines['left'].set_linewidth(bw)
ax0.spines['right'].set_linewidth(bw)
ax0.xaxis.set_tick_params(width=tw)
ax0.yaxis.set_tick_params(width=tw)
labels = getp(ax0, 'xticklabels')
for label in labels:
    label.set_fontsize(fontsize)
    label.set_fontname(fontname)
labels[0].set_visible(False)
labels[2].set_visible(False)
labels[4].set_visible(False)
labels[6].set_visible(False)
labels[8].set_visible(False)
labels[10].set_visible(False)
labels = getp(ax0, 'yticklabels')
for label in labels:
    label.set_fontsize(fontsize)
    label.set_fontname(fontname)

###############
# SECOND PLOT #
###############
ax1 = fig.add_axes(rect2, sharey=ax0)
for label in ax1.get_yticklabels():
    label.set_visible(False)
ax1.axis([0.05,0.55,epsmin,epsmax])

# draw plots
indicesR = (np.absolute(tempR-center_temp) < tol) & (np.absolute(np.log(rhoR/center_rho)) < tol)
indicesP = (np.absolute(tempP-center_temp) < tol) & (np.absolute(np.log(rhoP/center_rho)) < tol)
yetmp = yeP[indicesP]
epstmp = epsP[indicesP]
p_Pye, = ax1.plot(yetmp[1:],epstmp[1:],linecolor,linewidth=lw)
p_Rye, = ax1.plot(yeR[indicesR],epsR[indicesR],pointcolor,markersize=markersize,mew=mew)

ax1.spines['bottom'].set_linewidth(bw)
ax1.spines['top'].set_linewidth(bw)
ax1.spines['left'].set_linewidth(bw)
ax1.spines['right'].set_linewidth(bw)
ax1.xaxis.set_tick_params(width=tw)
ax1.yaxis.set_tick_params(width=tw)
labels = getp(ax1, 'xticklabels')
for label in labels:
    label.set_fontsize(fontsize)
    label.set_fontname(fontname)
labels = ax1.get_xticklabels()
for label in labels:
    label.set_fontsize(fontsize)
    label.set_fontname(fontname)
labels[0].set_visible(False)
minorLocator = MultipleLocator(0.02)
ax1.xaxis.set_minor_locator(minorLocator)

###############
# THIRD  PLOT #
###############
ax2 = fig.add_axes(rect1, sharey=ax0)
for label in ax2.get_yticklabels():
    label.set_visible(False)
ax2.axis([0.1,100,epsmin,epsmax])
ax2.set_xscale('log')

# draw plots
indicesR = (np.absolute(yeR-center_ye) < tol) & (np.absolute(np.log(rhoR/center_rho)) < tol)
indicesP = (np.absolute(yeP-center_ye) < tol) & (np.absolute(np.log(rhoP/center_rho)) < tol)
temptmp = tempP[indicesP]
epstmp = epsP[indicesP]
p_Ptemp, = ax2.plot(temptmp[1:],epstmp[1:],linecolor,linewidth=lw)
p_Rtemp, = ax2.plot(tempR[indicesR],epsR[indicesR],pointcolor,markersize=markersize,mew=mew)

ax2.spines['bottom'].set_linewidth(bw)
ax2.spines['top'].set_linewidth(bw)
ax2.spines['left'].set_linewidth(bw)
ax2.spines['right'].set_linewidth(bw)
ax2.xaxis.set_tick_params(width=tw)
ax2.yaxis.set_tick_params(width=tw)
labels = getp(ax2, 'xticklabels')
for label in labels:
    label.set_fontsize(fontsize)
    label.set_fontname(fontname)
labels = ax2.get_xticklabels()
for label in labels:
    label.set_fontsize(fontsize)
    label.set_fontname(fontname)
labels[1].set_visible(False)

## For some reason, this has to go after the plots
xlabelx = 0.5
xlabely = -0.2
ylabelx = -0.25
ylabely = 0.5

x0labelstring = r'$\rho\,[$g cm$^{-3}]$'
plt.text(xlabelx,xlabely,x0labelstring,fontsize=fontsize,
              horizontalalignment="center",transform=ax0.transAxes)

x1labelstring = "$Y_e$"
plt.text(xlabelx,xlabely,x1labelstring,fontsize=fontsize,
              horizontalalignment="center",transform=ax1.transAxes)

x2labelstring = r'$T\,[$MeV$]$'
plt.text(xlabelx,xlabely,x2labelstring,fontsize=fontsize,
              horizontalalignment="center",transform=ax2.transAxes)

Rhostring = r'$\rho=10^{10}\,\mathrm{g\,cm}^{-3}$'
Tstring   = r'$T=$'+str(center_temp)+' MeV'
Yestring  = r'$Y_e=$'+str(center_ye)
plt.text(0.13,0.85,Tstring  ,fontdict=font,transform=ax0.transAxes)
plt.text(0.1,0.72,Yestring ,fontdict=font,transform=ax0.transAxes)
plt.text(0.13,0.85,Tstring  ,fontdict=font,transform=ax1.transAxes)
plt.text(0.1345,0.72,Rhostring,fontdict=font,transform=ax1.transAxes)
plt.text(0.1,0.85,Yestring ,fontdict=font,transform=ax2.transAxes)
plt.text(0.1345,0.72,Rhostring,fontdict=font,transform=ax2.transAxes)




ylabelstring = r'$E\,[\mathrm{erg\,cm}^{-3}]$'
plt.text(ylabelx,ylabely,ylabelstring,fontdict=font,
         verticalalignment="center",rotation='vertical',transform=ax0.transAxes)

# legend
rc('legend', fontsize=fontsize)
ax0.legend([p_Prho,p_Rrho],
           ['Predicted','Solved'],
           loc=(0.05,0.35), handletextpad=0.1,prop=font).draw_frame(0)

plt.savefig('blackbody.pdf',bbox_inches='tight')#,transparent=True)


