import numpy as np
import math
import os
import sys
import glob
import shutil
import h5py

dataset_list = [
    #"annihilation_4force(erg|ccm|s,tet)",
    "four-force[abs](erg|ccm|s,tet)",
    "four-force[emit](erg|ccm|s,tet)",
    "l_abs(1|s|ccm,tet)",
    "l_emit(1|s|ccm,tet)",
    "distribution0(erg|ccm,tet)",
    "distribution1(erg|ccm,tet)",
    "distribution2(erg|ccm,tet)",
    "spectrum0(erg|s)",
    "spectrum1(erg|s)",
    "spectrum2(erg|s)"
]
print(dataset_list)

fluid_list = sorted(glob.glob("fluid_[0-9][0-9][0-9][0-9][0-9].h5"))[1:]
print(fluid_list)

# use file 0 as template
outfilename = "fluid_combined.h5"
#if not os.path.exists(outfilename):
shutil.copyfile("fluid_00001.h5",outfilename)
outfile = h5py.File(outfilename,"r+")

for dset in dataset_list:
    print(dset)

    # prepare datasets to combine
    data = np.array(outfile[dset])

    # loop through the fluid files
    for filename in fluid_list:
        infile = h5py.File(filename,"r")
        data += np.array(infile[dset])
        infile.close()

    outfile[dset][:] = data / len(fluid_list)

