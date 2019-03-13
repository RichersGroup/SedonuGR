import h5py
import numpy as np
import sys
sys.path.insert(0, '../')
import tools


rlname = "reflevel=0"

rho = 2e14 / tools.cactus_density # code units
temp = 10 # MeV
Ye = 0
gxx = 1
gyy = 2
gzz = 3
gxy = 0
gxz = 0
gyz = 0
lapse = 4


# NO SYMMETRIES
nx = 3
xmin = np.array([-1,-1,-1]) # Msun
xmax = np.array([ 1, 1, 1]) # Msun
delta = (xmax-xmin)/(nx-1)
extent = [xmin[0], xmax[0],
          xmin[1], xmax[1],
          xmin[2], xmax[2]]

outfile = h5py.File("stationary.h5","w")
outgroup = outfile.create_group(rlname)
outgroup.attrs["extent"] = extent
outgroup.attrs["delta"] = delta
outgroup.attrs["time"] = 0
outgroup.attrs["iteration"] = 0
outgroup.attrs["reflevel"] = 0

shape = np.ones((nx,nx,nx))
outgroup["rho"]   = rho   * shape
outgroup["temp"]  = temp  * shape
outgroup["Ye"]    = Ye    * shape
outgroup["gxx"]   = gxx   * shape
outgroup["gyy"]   = gyy   * shape
outgroup["gzz"]   = gzz   * shape
outgroup["gxy"]   = gxy   * shape
outgroup["gxz"]   = gxz   * shape
outgroup["gyz"]   = gyz   * shape
outgroup["lapse"] = lapse * shape
outgroup["velx"]  = 0     * shape
outgroup["vely"]  = 0     * shape
outgroup["velz"]  = 0     * shape
outgroup["betax"] = -np.array(outgroup["velx"])
outgroup["betay"] = -np.array(outgroup["vely"])
outgroup["betaz"] = -np.array(outgroup["velz"])

# SHIFTED
outfile_shifted = h5py.File("shifted.h5","w")
outfile.copy(rlname,outfile_shifted)
outgroup = outfile_shifted[rlname]
outgroup["velx"][...]  = .1 * shape
outgroup["vely"][...]  = .3 * shape
outgroup["velz"][...]  = .5 * shape
outgroup["betax"][...]  = -np.array(outgroup["velx"])
outgroup["betay"][...]  = -np.array(outgroup["vely"])
outgroup["betaz"][...]  = -np.array(outgroup["velz"])

outfile.close()
outfile_shifted.close()
