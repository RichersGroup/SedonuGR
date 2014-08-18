nx = 1
R_out = 1e5
dx = R_out/nx
rho = 1e-15
temp = 1e4;

r_in = 0
t0 = 100*3600.0*24
print '1D_sphere',nx,r_in

for i  in range(1,nx+1):
    print i*dx,rho,temp,
    print 1.0


