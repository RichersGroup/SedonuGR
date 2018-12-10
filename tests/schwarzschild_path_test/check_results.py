import h5py
import numpy as np
import matplotlib.pyplot as plt

tolerance = .05
k_b = 1.3806488e-16 # erg/K
c = 2.99e10         #cm/s
h = 6.62606957e-27  #erg/Hz
pi = 3.14159265359
k_MeV = 1.16046e10
MeV = 0.0000016021773
r = 1.5 #cm
T = 10*k_MeV #K
mue = 10*MeV
rs = 1.0
alpha_core = np.sqrt(1.-rs/r)
Rout = 10

kr_0 = 1
kt_0 = 1

# do error checking
with open("radial.dat","r") as search:
    for line in search:
        if "DO_GR" in line:
            do_gr = int(line[-2])

radial_k = []
radial_k.append(np.genfromtxt("radial.dat",usecols=(4))[-1])
radial_k.append(np.genfromtxt("radial.dat",usecols=(5))[-1])
radial_k.append(np.genfromtxt("radial.dat",usecols=(6))[-1])
radial_k.append(np.genfromtxt("radial.dat",usecols=(7))[-1])
radial_x = []
radial_x.append(np.genfromtxt("radial.dat",usecols=(0))[-1])
radial_x.append(np.genfromtxt("radial.dat",usecols=(1))[-1])
radial_x.append(np.genfromtxt("radial.dat",usecols=(2))[-1])

around_k = []
around_k.append(np.genfromtxt("around.dat",usecols=(4))[-1])
around_k.append(np.genfromtxt("around.dat",usecols=(5))[-1])
around_k.append(np.genfromtxt("around.dat",usecols=(6))[-1])
around_k.append(np.genfromtxt("around.dat",usecols=(7))[-1])
around_x = []
around_x.append(np.genfromtxt("around.dat",usecols=(0))[-1])
around_x.append(np.genfromtxt("around.dat",usecols=(1))[-1])
around_x.append(np.genfromtxt("around.dat",usecols=(2))[-1])

if(do_gr==0):
    radial_k_expected = [1,0,0,1]
    radial_x_expected = [Rout,0,0]
    around_k_expected = [0,1,0,1]
    around_x_expected = [r,np.sqrt(Rout**2-r**2),0]
else:
    radial_k_expected = [1-rs/r,0,0, (1.-rs/r) / (1.-rs/radial_x[0]) ]
    radial_x_expected = [Rout,0,0]
    around_k_expected = [-np.sqrt(1-rs/r),0,0,1]
    around_x_expected = [0,r,0]

print(radial_k_expected, radial_k)
print(radial_x_expected, radial_x)
print(around_k_expected, around_k)
print(around_x_expected, around_x)
    
for i in range(4):
    radial_error = abs(radial_k_expected[i] - radial_k[i])
    around_error = abs(around_k_expected[i] - around_k[i])
    if radial_error>tolerance or around_error>tolerance:
        raise Exception("spherical_emis results are outside of the tolerance.")

for i in range(3):
    radial_error = abs(radial_x_expected[i] - radial_x[i])
    around_error = abs(around_x_expected[i] - around_x[i])
    if radial_error>tolerance or around_error>tolerance:
        raise Exception("spherical_emis results are outside of the tolerance.")

print("SUCCESS")
exit()

