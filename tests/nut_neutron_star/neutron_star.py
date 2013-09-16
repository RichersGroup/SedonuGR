MeV = 1.1605e10 # Kelvin
km = 1.0e5      # cm

nx = 100
t0 = -1

R_max = 12*km
R_min = 11.5*km
dx = (R_max-R_min)/nx
rho_max = 2e15
ye_max = 0.2

temp = 5*MeV

print '1D_sphere', 'GRB',nx,R_min,t0

for i  in range(1,nx+1):
    R = R_min + i*dx
    print R, rho_max*(1-(R/R_max)**2), temp, ye_max*(1-(R/R_max)**2)+.05
