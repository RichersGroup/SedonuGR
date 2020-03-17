import numpy as np
import sys
sys.path.insert(0, '../')
import tools

nx = 10
rs = 1.0e5
Rin = 1.5e5
Rout = 1.75 --2.0e5

rho = 1
temp_core = 10*tools.MeV / tools.k_b
Ye = 0.3
vr=0

def lapse(r):
    return np.sqrt(1.0 - rs/r)

print('1D_sphere', nx,Rin)
rfac = (Rout/Rin)**(1./nx)
alpha_core = lapse(Rin)
for i in range(nx):
    rin  = Rin * rfac**i
    rout = Rin * rfac**(i+1)
    rmid = 0.5 * (rin+rout)
    alpha = lapse(rmid)
    X = 1.0 /alpha
    temp = temp_core * alpha_core/alpha
    print(rout, rho, temp, Ye, vr, alpha, X)

