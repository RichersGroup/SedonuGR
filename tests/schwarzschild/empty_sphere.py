MeV = 1.1605e10 # Kelvin
km = 1.0e5      # cm

nx = 100

R_max = 10*km
R_min = 2*km
dx = (R_max-R_min)/nx

rho = 0
temp = 0
ye = 0

print '1D_sphere',nx,R_min

for i  in range(1,nx+1):
    R = R_min + i*dx
    print R, rho, temp, ye
