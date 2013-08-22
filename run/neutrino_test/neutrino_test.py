MeV = 1.1605e10
km = 1.0e5

nx = 100
t0 = -1
v_in = 0

dx = 12*km/nx
temp = 30*MeV
rho = 3e14

print '1D_sphere', 'GRB',nx,v_in,t0

for i  in range(1,nx+1):
    print i*dx,rho,temp,
    print 1.0


