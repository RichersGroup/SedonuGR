import numpy as np

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

#kup_initial = [1,1,0,0] # traveling in the x-direction
mean_kup_expected = np.array([0,0,0,1]) * 8.01089e-07

# uncomment data reading when data is output.
kup_data = np.genfromtxt("elastic_isotropic_kernel.dat") #[particle, component]
#kup_data = 2*np.random.rand(10000,4)-1
#kup_data[:,0] = 1

nparticles = np.shape(kup_data)[0]
kup_mean = np.sum(kup_data,axis=0)/nparticles
kup_mean2 = np.sum(kup_data*kup_data, axis=0)/nparticles
kup_stddev = np.sqrt(np.abs(kup_mean2 - kup_mean**2))/nparticles**.25

# check that mean is within standard deviation of expected value
for i in range(3):
    if(np.abs(kup_mean[i] - mean_kup_expected[i]) > kup_stddev[i]):
        raise Exception("Expected:"+str(mean_kup_expected)+" but got "+str(kup_mean)+"+/-"+str(kup_stddev)+" (component "+str(i)+")")

print("SUCCESS: Expected "+str(mean_kup_expected)+" and got "+str(kup_mean)+"+/-"+str(kup_stddev))
exit()

