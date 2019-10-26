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
gyy = 1 #2
gzz = 1 #3
gxy = 0
gxz = 0
gyz = 0
lapse = 1 #4


# NO SYMMETRIES
nx = 4
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
outgroup["velx"]  =  0.2     * shape
outgroup["vely"]  = -0.4     * shape
outgroup["velz"]  =  0.6     * shape
outgroup["betax"] = -0.1 * shape
outgroup["betay"] =  0.3 * shape
outgroup["betaz"] = -0.5 * shape

outfile.close()
