import matplotlib.pyplot as plt
import h5py
import numpy as np
import os

############
# PLOTTING #
############
infilename = "0.333_s012_g0-16_nf9_LxxxLyyyLzzz"
f = h5py.File("analysis/"+infilename+".h5","r")
xgrid = np.array(f["xgrid"])
ygrid = np.array(f["ygrid"])
zgrid = np.array(f["zgrid"])

outdir = "plots/"+infilename
if not os.path.exists(outdir):
    os.makedirs(outdir)

for normalDirection in range(3):
    print("direction",normalDirection)
    if normalDirection==0:
        xlabel = "z (km)"
        ylabel = "y (km)"
        plotdelta_2D = np.array(f["err_yz"])
        plotX2_2D    = np.array(f["chi2_yz"])
        grid1 = zgrid
        grid2 = ygrid
    if normalDirection==1:
        xlabel = "z (km)"
        ylabel = "x (km)"
        plotdelta_2D = np.array(f["err_xz"])
        plotX2_2D    = np.array(f["chi2_xz"])
        grid1 = zgrid
        grid2 = xgrid
    if normalDirection==2:
        xlabel = "y (km)"
        ylabel = "x (km)"
        plotdelta_2D = np.array(f["err_xy"])
        plotX2_2D    = np.array(f["chi2_xy"])
        grid1 = ygrid
        grid2 = xgrid

    plotdelta_2D[np.where(plotdelta_2D!=plotdelta_2D)] = -1.0
    plotX2_2D[   np.where(   plotX2_2D!=   plotX2_2D)] = -1.0

    plt.clf()
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.pcolormesh(grid1/1e5,grid2/1e5,plotdelta_2D,vmin=0)
    plt.axes().set_aspect('equal', 'datalim')
    plt.colorbar(label=r"$\delta$")
    plt.axes().axis('tight')
    plt.savefig(outdir+"/delta_"+str(normalDirection)+".png")

    plt.clf()
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.axes().set_aspect('equal', 'datalim')
    plt.pcolormesh(grid1/1e5,grid2/1e5,plotX2_2D,vmin=0)
    plt.colorbar(label=r"$\chi^2/\nu$")
    plt.axes().axis('tight')
    plt.savefig(outdir+"/X2_"+str(normalDirection)+".png")
