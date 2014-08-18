nx = 100
vmax = 1e9
t0 = 100*3600.0*24
dx = vmax*t0/nx
r_in = 0
xmax = vmax*t0
temp = 1e4;

rho = 1.4*1.99e33/(4.0*3.14/3.0*xmax**3.0)

print '1D_sphere',nx,r_in

for i  in range(1,nx+1):
    print i*dx,rho,temp,
    print 1.0


