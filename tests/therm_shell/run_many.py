import os
import subprocess
from math import *
from scipy.integrate import quad
from numpy import Inf

MeV = 1.1605e10 # Kelvin
km = 1.0e5      # cm

# INPUTS
eosfile   = "../../external/tables/EOS/SFHo.h5"

R_min = 10*km
R_frac = 1.e-4

min_logrho = 6  #g/ccm
max_logrho = 15 #g/ccm
n_rho = 10
center_logrho = 10

min_logT = -1 #MeV
max_logT = 2  #MeV
n_T = 10
center_logT = 1

min_ye = 0.05
max_ye = 0.55
n_ye = 10
center_ye = 0.3


# create the shell
nx = 1

R_max = R_min*(1.+R_frac)
dx = (R_max-R_min)/nx

rho = 1e10
temp = 10*MeV
ye = 0.3
alpha = 1
X = 1

f = open("shell.mod","w")
print('1D_sphere',nx,R_min,file=f)
print(R_max, rho, temp, ye,0,alpha,X, file=f)
f.close()


# set parameter file
string = "sed " + \
         "-e 's/RCORE_HERE/"    + '%.5e' % R_min         + "/g' " + \
         "param.template > param.lua"
os.system(string)

def plot(rho0,T0,ye0):
    string = "sed " + \
             "-e 's/RHO_HERE/"    + '%.5e' % rho0         + "/g' " + \
             "-e 's/TEMP_HERE/"   + str(T0)           + "/g' " + \
             "-e 's/YE_HERE/"     + str(ye0)          + "/g' " + \
             "template.gnuplot > compare.gnuplot"
    os.system(string)
    os.system("gnuplot compare.gnuplot")


### BEGIN SCRIPT ###

rho0 = 10**(center_logrho)
T0   = 10**(center_logT)
ye0  =      center_ye
dlogrho = (max_logrho - min_logrho) / (n_rho - 1.0)
dlogT   = (max_logT   - min_logT  ) / (n_T   - 1.0)
dye     = (max_ye     - min_ye    ) / (n_ye  - 1.0)

string = "../../exe/therm_shell_test param.lua " + \
         str(min_logrho) + " " + str(max_logrho) + " " + str(rho0) + " " + str(n_rho) + " " + \
         str(min_logT  ) + " " + str(max_logT  ) + " " + str(T0  ) + " " + str(n_T  ) + " " + \
         str(min_ye    ) + " " + str(max_ye    ) + " " + str(ye0 ) + " " + str(n_ye ) + " " + \
         eosfile
print(string)
os.system(string)

plot(rho0,T0,ye0)
