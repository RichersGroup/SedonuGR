#!/usr/bin/python                                                               
import sys

from numpy import *
import numpy as np
from astropy.io import ascii

if len(sys.argv)!=3:
    print("Usage: python short_to_mod.py file.short interpolation_multiplier")
    sys.exit()

data = loadtxt(sys.argv[1],usecols=[2,3,4,6],skiprows=1)
interpolation_multiplier = int(sys.argv[2])
r_f   = data[:,0]
T_f   = data[:,1]
rho_f = data[:,2]
Ye_f  = data[:,3]
print len(r_f)
r   = np.ones( (len(r_f)-1) * interpolation_multiplier + 1 )
T   = np.ones( (len(r_f)-1) * interpolation_multiplier + 1 )
rho = np.ones( (len(r_f)-1) * interpolation_multiplier + 1 )
Ye  = np.ones( (len(r_f)-1) * interpolation_multiplier + 1 )
print len(r)

# first zone is unchanged, constant rho,T,Ye
r[0]   =   r_f[0]
rho[0] = rho_f[0]
Ye[0]  =  Ye_f[0]
T[0]   =   T_f[0]

# logarithmic interpolation
rL = 0
rR = 0
for i in range(1,len(r_f)):
    r_ratio   = pow(   r_f[i]/  r_f[i-1], 1.0/interpolation_multiplier)
    T_ratio   = pow(   T_f[i]/  T_f[i-1], 1.0/interpolation_multiplier)
    Ye_ratio  = pow(  Ye_f[i]/ Ye_f[i-1], 1.0/interpolation_multiplier)
    rho_ratio = pow( rho_f[i]/rho_f[i-1], 1.0/interpolation_multiplier)

    rho_power = log(rho_ratio)/log(r_ratio)
    Ye_power  = log( Ye_ratio)/log(r_ratio)
    T_power   = log(  T_ratio)/log(r_ratio)
    eps_power = T_power + rho_power

    for j in range(0,interpolation_multiplier):
        index = 1 + j + (i-1)*interpolation_multiplier
        rL   =   r_f[i-1] * pow(  r_ratio,j)
        TL   =   T_f[i-1] * pow(  T_ratio,j)
        YeL  =  Ye_f[i-1] * pow( Ye_ratio,j)
        rhoL = rho_f[i-1] * pow(rho_ratio,j)

        # ratios are defined assuming conservative variables (rho,eps,Ye) scale as some power of r. Assume eps = nkT.
        # integrated by hand to average over cell volume
        r[index]   = r_f[i-1] * pow(r_ratio, j+1)
        Ye[index]  = 3* YeL/(3+ Ye_power) * (pow(r_ratio, Ye_power+3)-1) / (pow(r_ratio,3)-1)
        rho[index] = 3*rhoL/(3+rho_power) * (pow(r_ratio,rho_power+3)-1) / (pow(r_ratio,3)-1)
        T[index]   = 3*  TL/(3+eps_power) * (pow(r_ratio,eps_power+3)-1) / (pow(r_ratio,3)-1) * rhoL/rho[index]

savetxt(sys.argv[1]+".mod", array([r,rho,T,Ye]).T, header="1D_sphere %d 0"%len(r), comments='')
