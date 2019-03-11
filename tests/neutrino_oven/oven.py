import numpy as np
import sys
sys.path.insert(0, '../')
import tools


nx = 10
rs = 1.0e5
Rin = 1.5e5
Rout = 2.0e5

rho = 1
temp_core = 10*tools.MeV / tools.k_b
Ye = 0.3
vr=0

print('1D_sphere', nx,Rin)
rfac = (Rout/Rin)**(1./nx)
alpha_core = np.sqrt(1.0 - rs/Rin)
for i in range(nx):
    rin  = Rin * rfac**i
    rout = Rin * rfac**(i+1)
    rmid = 0.5 * (rin+rout)
    alpha = alpha_core #np.sqrt(1.0-rs/rmid)
    X = 1.0 #/alpha
    temp = temp_core * (alpha/alpha_core)**0 #(alpha_core/alpha)
    print(rout, rho, temp, Ye, vr, alpha, X)

