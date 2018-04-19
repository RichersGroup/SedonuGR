import numpy as np
nx = 9*4
rs = 1.0
Rin = 1
Rout = 10

rho = 0
temp = 0
Ye = 0
vr=0

print('1D_sphere', nx,Rin)
rfac = (Rout/Rin)**(1./nx)
for i in range(nx):
    rin  = Rin * rfac**i
    rout = Rin * rfac**(i+1)
    rmid = 0.5 * (rin+rout)
    alpha = np.sqrt(1.0-rs/rmid)
    X = 1.0/alpha
    print(rout, rho, temp, Ye, vr, alpha, X)

