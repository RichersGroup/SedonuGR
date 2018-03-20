import numpy as np
import math
import os
import sys
import glob
import shutil
import h5py

# use file 0 as template
outfilename = "fluid_combined.h5"
#if not os.path.exists(outfilename):
shutil.copyfile("fluid_00000.h5",outfilename)
outfile = h5py.File(outfilename,"r+")

# prepare datasets to combine
force_annihil = np.zeros(np.shape(outfile["annihilation_4force(erg|ccm|s,tet)"]))
force_abs     = np.zeros(np.shape(outfile["four-force[abs](erg|ccm|s,tet)"]))
force_emit    = np.zeros(np.shape(outfile["four-force[emit](erg|ccm|s,tet)"]))
l_abs         = np.zeros(np.shape(outfile["l_abs(1|s|ccm,tet)"]))
l_emit        = np.zeros(np.shape(outfile["l_emit(1|s|ccm,tet)"]))
distribution0 = np.zeros(np.shape(outfile["distribution0(erg|ccm,tet)"]))
distribution1 = np.zeros(np.shape(outfile["distribution1(erg|ccm,tet)"]))
distribution2 = np.zeros(np.shape(outfile["distribution2(erg|ccm,tet)"]))
spectrum0     = np.zeros(np.shape(outfile["spectrum0(erg|s)"]))
spectrum1     = np.zeros(np.shape(outfile["spectrum1(erg|s)"]))
spectrum2     = np.zeros(np.shape(outfile["spectrum2(erg|s)"]))

# loop through the fluid files
fluid_list = sorted(glob.glob("fluid_[0-9][0-9][0-9][0-9][0-9].h5"))[1:]
print(fluid_list)
for filename in fluid_list:
    infile = h5py.File(filename,"r")
    force_annihil += np.array(infile["annihilation_4force(erg|ccm|s,tet)"])
    force_abs     += np.array(infile["four-force[abs](erg|ccm|s,tet)"])
    force_emit    += np.array(infile["four-force[emit](erg|ccm|s,tet)"])
    l_abs         = np.array(infile["l_abs(1|s|ccm,tet)"])
    l_emit        += np.array(infile["l_emit(1|s|ccm,tet)"])
    distribution0 += np.array(infile["distribution0(erg|ccm,tet)"])
    distribution1 += np.array(infile["distribution1(erg|ccm,tet)"])
    distribution2 += np.array(infile["distribution2(erg|ccm,tet)"])
    spectrum0     += np.array(infile["spectrum0(erg|s)"])
    spectrum1     += np.array(infile["spectrum1(erg|s)"])
    spectrum2     += np.array(infile["spectrum2(erg|s)"])
    print(filename)
print(np.sum(l_abs))
# write out normalized files
outfile["annihilation_4force(erg|ccm|s,tet)"][:] = force_annihil / len(fluid_list)
outfile["four-force[abs](erg|ccm|s,tet)"][:]     = force_abs     / len(fluid_list)
outfile["four-force[emit](erg|ccm|s,tet)"][:]    = force_emit    / len(fluid_list)
outfile["l_abs(1|s|ccm,tet)"][:]                 = l_abs         / len(fluid_list)
outfile["l_emit(1|s|ccm,tet)"][:]                = l_emit        / len(fluid_list)
outfile["distribution0(erg|ccm,tet)"][:]         = distribution0 / len(fluid_list)
outfile["distribution1(erg|ccm,tet)"][:]         = distribution1 / len(fluid_list)
outfile["distribution2(erg|ccm,tet)"][:]         = distribution2 / len(fluid_list)
outfile["spectrum0(erg|s)"][:]                   = spectrum0     / len(fluid_list)
outfile["spectrum1(erg|s)"][:]                   = spectrum1     / len(fluid_list)
outfile["spectrum2(erg|s)"][:]                   = spectrum2     / len(fluid_list)
