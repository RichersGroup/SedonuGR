MeV = 1.1605e10
km = 1.0e5

nx = 100
t0 = -1
v_in = 0

R_max = 10*km
R_min = 9*km
dx = (R_max-R_min)/nx
rho_max = 2e15
ye_max = 0.2

temp = 1*MeV

print '1D_sphere', 'GRB',nx,v_in,t0

for i  in range(1,nx+1):
    R = R_min + i*dx
    print R, rho_max*(1-(R/R_max)**2), temp, ye_max*(1-(R/R_max)**2)+.05
