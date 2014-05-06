MeV = 1.1605e10 # Kelvin
km = 1.0e5      # cm

nx = 1
t0 = -1

R_max = 10000*km
R_min = 0*km
dx = (R_max-R_min)/nx

rho = RHOHERE
temp = TEMPHERE*MeV
ye = YEHERE

print '1D_sphere', 'GRB',nx,R_min,t0

for i  in range(1,nx+1):
    R = R_min + i*dx
    print R, rho, temp, ye
