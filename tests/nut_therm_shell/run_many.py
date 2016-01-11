import os
import subprocess
from math import *
from scipy.integrate import quad
from numpy import Inf

# INPUTS
eosfile   = "/data/tables/EOS/HShen.h5"

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

def plot(rho0,T0,ye0):
    string = "sed " + \
             "-e 's/RHO_HERE/"    + str(rho0)         + "/g' " + \
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

os.system("python shell.py > shell.mod")
string = "mpirun -np 2 -env OMP_NUM_THREADS 2 ./nut_therm_shell.exe param.lua " + \
         str(min_logrho) + " " + str(max_logrho) + " " + str(rho0) + " " + str(n_rho) + " " + \
         str(min_logT  ) + " " + str(max_logT  ) + " " + str(T0  ) + " " + str(n_T  ) + " " + \
         str(min_ye    ) + " " + str(max_ye    ) + " " + str(ye0 ) + " " + str(n_ye ) + " " + \
         eosfile
print string
os.system(string)

plot(rho0,T0,ye0)
