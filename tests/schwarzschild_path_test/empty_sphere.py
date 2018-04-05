import numpy as np
nx = 9*8
rs = 1.0
rin = 1

R = 10
rho = 0
temp = 0
Ye = 0
vr=0

print('1D_sphere', nx,rin)
dr = (R-rin)/nx
for i in range(nx):
    rout = rin + dr*(i+1  )
    rmid = rin + dr*(i+0.5)
    alpha = np.sqrt(1.0-rs/rmid)
    X = 1.0/alpha
    print(rout, rho, temp, Ye, vr, alpha, X)

