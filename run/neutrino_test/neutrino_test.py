MeV = 1.1605e10
km = 1.0e5

nx = 100
t0 = -1
v_in = 0

R = 12*km
rho_max = 2e12
dx = R/nx

temp = 10*MeV

print '1D_sphere', 'GRB',nx,v_in,t0

for i  in range(1,nx+1):
    print i*dx, rho_max*(1-(float(i)/float(nx))**2), temp, 1.0
