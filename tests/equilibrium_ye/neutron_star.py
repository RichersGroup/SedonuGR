MeV = 1.1605e10
km = 1.0e5

nx = 100

R_max = 10*km
R_min = 9*km
dx = (R_max-R_min)/nx
rho_max = 2e12
ye_max = 0.2

temp = 5*MeV
lapse = 1.
X = 1.

print('1D_sphere',nx,R_min)

for i  in range(1,nx+1):
    R = R_min + i*dx
    print(R, rho_max*(1-(R/R_max)**2), temp, ye_max,0,lapse,X)
