import numpy as np
import sys
sys.path.insert(0, '../')
import tools


nx = 10
rs = 1.0e5
Rin = 1.5e5
Rout = 2.0e5

rho = 1
temp = 10*tools.MeV / tools.k_b
Ye = 0.3
vr=0

#with open("output.txt") as search:
#    for line in search:
#        if "DO_GR" in line:
#            if(int(line[-2])==0):
#                alpha = np.ones(len(r))
 #               alpha_core = 1
#            else:
#                alpha = np.sqrt(1. - rs/r)
#                alpha_core = np.sqrt(1.-rs/r0)
#T_theory = T0/k_MeV * (alpha_core/alpha)**3


print('1D_sphere', nx,Rin)
rfac = (Rout/Rin)**(1./nx)
for i in range(nx):
    rin  = Rin * rfac**i
    rout = Rin * rfac**(i+1)
    rmid = 0.5 * (rin+rout)
    alpha = np.sqrt(1.0-rs/rmid)
    X = 1.0/alpha
    print(rout, rho, temp, Ye, vr, alpha, X)

