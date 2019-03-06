import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import math
import h5py

sigma = 200

# read in data
filename = "fluid_00001.h5"
f = h5py.File(filename, "r")
edens = np.array(f.get("distribution0(erg|ccm,tet)")[:,0,0])[1:] # r, e, moment
rgrid = f.get("axes/x0(cm)[edge]")[1:]
rcenter = f.get("axes/x0(cm)[mid]")[1:]
flux  = np.array(f.get("distribution0(erg|ccm,tet)")[:,0,1])[1:] # r, e, moment
nr = len(rcenter)

## binomial smooth edens
#edens = np.zeros(np.shape(edens_unsmoothed))
#edens[0] = (2.*edens_unsmoothed[0] + edens_unsmoothed[1])/3.
#edens[-1] = (edens_unsmoothed[-2] + 2.*edens_unsmoothed[-1])/3.
#for i in range(1,len(edens)-1):
#    edens[i] = (edens_unsmoothed[i-1] + 2.*edens_unsmoothed[i] + edens_unsmoothed[i+1])/4.

# set up the grid
rface = rgrid[1:-1]

# evaluate Fick's law
gradE = np.zeros(nr-1)
gradE[0] = (edens[1]-edens[0]) / (rcenter[1]-rcenter[0])
gradE[-1] = (edens[-1]-edens[-2]) / (rcenter[-1]-rcenter[-2])
for i in range(1, nr-1):
    gradE[i] = (edens[i+1]-edens[i-1]) / (rcenter[i+1]-rcenter[i-1])
fick_flux = -gradE / (3.0 * sigma)

# get chi2
flux_face = np.array([ (flux[i]+flux[i+1])/2. for i in range(len(fick_flux))])
mean = np.sum(fick_flux)/len(fick_flux)
mean2 = np.sum(fick_flux**2)/len(fick_flux)
stddev = np.sqrt(mean2 - mean**2)
error2 = (flux_face-fick_flux)**2/stddev**2
chi2 = np.sum(error2) / len(flux)
print(mean, stddev, chi2)

meanflux = np.sum(flux)/len(flux)
if(np.abs(mean-meanflux)/stddev > 1):
    raise Exception("Mean_flux="+str(meanflux)+" does not agree with mean_fick_flux="+str(mean)+" within the stddev=+/-"+str(stddev))
else:
    print("SUCCESS")
    print("Mean_flux="+str(meanflux)+" agrees with mean_fick_flux="+str(mean)+" within the stddev=+/-"+str(stddev))


## Set up plots
panelheight = 0.98  # height of a single panel
panelheighttop = 0.7  # height of a single panel
left, bottom, width, height = 0.15, 0.15, 0.8, 0.8
rect0 = [left, bottom, width, height*panelheight]
fig = plt.figure()#figsize=(11., 12.))
font = 20

# annotation text
fig.gca().get_xaxis().set_visible(False)
fig.gca().get_yaxis().set_visible(False)
plt.setp(fig.gca(),frame_on=False)

min1 = min(flux_face)
min2 = min(fick_flux)
max1 = max(flux_face)
max2 = max(fick_flux)
ax0 = fig.add_axes(rect0)
ax0.axis([min(rgrid), max(rgrid), min(min1,min2), max(max1,max2)])#0, 4e22])#

# draw plots
plist = []
p, = ax0.plot(rface,flux_face, color="k", linewidth=1)
plist.append(p)
p, = ax0.plot(rface,fick_flux, color="b", linewidth=1)
plist.append(p)
ax0.axhline(mean, color="gray")
ax0.axhline(mean+stddev, color="gray")
ax0.axhline(mean-stddev, color="gray")

ax0.legend(plist, ["$\mathrm{Flux}\,[\mathrm{erg/ccm}]$",r"$-\nabla E / 3\sigma\,[\mathrm{erg/ccm}]$"],
           loc=(0.5,0.8),).draw_frame(0)

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
