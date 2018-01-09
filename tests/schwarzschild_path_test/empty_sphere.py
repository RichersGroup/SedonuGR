import numpy as np
nx = 999
rs = 1.0
rin = 1.34

R = 10
rho = 0
temp = 0
Ye = 0
vr=0

print('1D_sphere', nx,rin)
for i in range(nx):
    rout = rin + (R-rin)/nx*(i+1  )
    rmid = rin + (R-rin)/nx*(i+0.5)
    alpha = np.sqrt(1.0-rs/rmid)
    X = 1.0/alpha
    print(rout, rho, temp, Ye, vr, alpha, X)

