import numpy as np
import sys

nx = 25
nx_15 = int(sys.argv[1])
rs = 1.0
Rin = 1

rho = 0
temp = 0
Ye = 0
vr=0

print('1D_sphere', nx,Rin)
rfac = (1.5)**(1./nx_15)
for i in range(nx):
    rin  = Rin * rfac**i
    rout = Rin * rfac**(i+1)
    rmid = 0.5 * (rin+rout)
    alpha = np.sqrt(1.0-rs/rmid)
    X = 1.0/alpha
    print(rout, rho, temp, Ye, vr, alpha, X)

