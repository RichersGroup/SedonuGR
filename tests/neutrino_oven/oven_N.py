import numpy as np
import sys
sys.path.insert(0, '../')
import tools

nx = 31
rs = 1.0e5
Rcore = 1.5e5
Rout = 3e5 #1.75 --2.0e5

rho = 1
temp_core = 10*tools.MeV / tools.k_b
Ye = 0.3
vr=0

def lapse(r):
    return np.sqrt(1.0 - rs/r)

print('1D_sphere', nx-2,0)
dr = Rout/nx
alpha_core = lapse(Rcore)
for i in range(2,nx):
    rin  = dr*i
    rout = dr*(i+1)
    rmid = 0.5 * (rin+rout)
    if(rout<=Rcore):
        alpha = alpha_core
    else:
        alpha = lapse(rmid)
    X = 1.0/alpha
    temp = temp_core #* alpha_core/alpha
    print(rout, rho, temp, Ye, vr, alpha, X)

