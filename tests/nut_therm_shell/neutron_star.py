MeV = 1.1605e10 # Kelvin
km = 1.0e5      # cm

nx = 1
t0 = -1

R_max = 10*km
R_min = 9.9999*km
dx = (R_max-R_min)/nx

rho = 0 #1e12
temp = 10*MeV
ye = 0.181845

print '1D_sphere', 'GRB',nx,R_min,t0

for i  in range(1,nx+1):
    R = R_min + i*dx
    print R, rho, temp, ye
