import numpy as np
import math
import os
import sys
import glob
import shutil
import h5py

# must have decimal so python interprets as float and not int
# face-centered input values
new_xmin = [0.,0.,0.] # Msun
new_xmax = [2.,2.,2.] # Msun

infilename = "/panfs/ds08/sxs/scurtis/Slices/SFHo_M140120_M0/938650.h5"
infile = h5py.File(infilename, "r")
reflevel = 4
rlname = "reflevel="+str(reflevel)

ingroup = infile[rlname]
delta = ingroup.attrs["delta"]
extent = ingroup.attrs["extent"]
xmin = [extent[0],extent[2],extent[4]]
xmax = [extent[1],extent[3],extent[5]]
shape = np.shape(np.array(ingroup["Ye"]))
print("delta: ",delta)
print("old extent: ",extent)

# check that inputs make sense
for i in range(3):
    assert(new_xmin[i] >= extent[2*i  ] - delta[i]/2.)
    assert(new_xmax[i] <= extent[2*i+1] + delta[i]/2.)

xvertex = [[xmin[i] + delta[i]*j for j in range(shape[i])] for i in range(len(shape))]
imin = [np.where(np.array(xvertex[i]-delta[i]/2.)<=new_xmin[i])[0][-1] for i in range(len(shape))]
imax = [np.where(np.array(xvertex[i]+delta[i]/2.)>=new_xmax[i])[0][0] for i in range(len(shape))]

outfilename = "domain_cut.h5"
outfile = h5py.File(outfilename,"w")
outgroup = outfile.create_group(rlname)

# new extent should be vertex-centered
new_extent = [
    xmin[0] + imin[0]*delta[0], xmin[0] + imax[0]*delta[0],
    xmin[1] + imin[1]*delta[1], xmin[1] + imax[1]*delta[1],
    xmin[2] + imin[2]*delta[2], xmin[2] + imax[2]*delta[2]]
print("new extent: ",new_extent)
print("requested domain: ",new_xmin, new_xmax)
print("actual domain: ",
      [xmin[i] + (imin[i]-.5)*delta[i] for i in range(3)], 
      [xmin[i] + (imax[i]+.5)*delta[i] for i in range(3)])
outgroup.attrs["extent"] = new_extent
copy_attributes = ["delta","time","iteration","reflevel"]
for attr in copy_attributes:
    outgroup.attrs[attr] = ingroup.attrs[attr]

dataset_list = [
    "Ye",
    "betax","betay","betaz",
    "entr",
    "eps",
    "gxx","gxy","gxz","gyy","gyz","gzz",
    "lapse",
    "press",
    "rho",
    "temp",
    "velx","vely","velz",
    "vol",
    "w_lorentz"]
for dset in dataset_list:
    outgroup[dset] = np.array(ingroup[dset])[imin[0]:imax[0]+1, imin[1]:imax[1]+1, imin[2]:imax[2]+1]

outfile.close()
infile.close()
