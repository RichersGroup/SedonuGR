import numpy as np
import sys
sys.path.insert(0, '../')
import tools

nx = 500
R = 1
Rout = 2

rho = 1
T = 10*tools.MeV / tools.k_b
Ye = 0.3
vr=0

print('1D_sphere', nx-2,0)
dr = Rout/nx
for i in range(2,nx):
    rin  = dr*i
    rout = dr*(i+1)
    rmid = 0.5 * (rin+rout)
    alpha = 1
    X = 1.0/alpha

    if(rmid < R):
        print(rout, rho, T, Ye, vr, alpha, X)
    else:
        print(rout, 0  , 0   , 0 , vr, alpha, X)

