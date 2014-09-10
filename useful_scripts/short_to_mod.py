#!/usr/bin/python                                                               
import sys

from numpy import *
from astropy.io import ascii

if len(sys.argv)!=2:
    print("Usage: python short_to_mod.py file.short")
    sys.exit()

data = loadtxt(sys.argv[1],usecols=[2,3,4,6],skiprows=1)

r   = data[:,0]
T   = data[:,1]
rho = data[:,2]
Ye  = data[:,3]

savetxt(sys.argv[1]+".mod", array([r,rho,T,Ye]).T, header="1D_sphere %d 0"%len(r), comments='')
