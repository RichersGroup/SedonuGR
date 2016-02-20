import os
from math import *
from scipy.integrate import quad
import numpy as np

# INPUTS
eosfile   = "/data/tables/EOS/HShen.h5" #ignored if compiled for helmholtz eos

min_logrho = 8  #g/ccm
max_logrho = 15 #g/ccm
n_rho = 100
center_logrho = 10

min_logT = -1 #MeV
max_logT = 2  #MeV
n_T = 100
center_logT = 0.5

min_ye = 0.05
max_ye = 0.55
n_ye = 100
center_ye = 0.3


# Constants
h = 6.626075518e-27 #erg s
c = 2.99792458e10   #cm/s
MeV = 1.60217657e-6  #erg/MeV


def BB(x,T,mu):
    # x = h nu / kT
    # T and mu are in same units
    return x**3 / ( exp(x - mu/T) + 1 )

def edens(T,mu):
    # mu and T are in the same units
    if(mu/T > 10): max_x = 3*mu/T
    else: max_x = 20
    return 4*pi*h*c * (T*MeV/c/h)**4 * quad(BB,0.0,max_x,args=(T,mu))[0]

def write_predicted(mu,rho,T,Ye):
    etot = edens(T,mu) + edens(T,-mu) + 4*edens(T,0)
    string = str(etot)+" "+str(rho)+" "+str(T)+" "+str(Ye)
    os.system("echo "+string+" >> predicted.dat")

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
T0   = 10**(center_logT  )
ye0  =      center_ye
dlogrho = (max_logrho - min_logrho) / (n_rho - 1.0)
dlogT   = (max_logT   - min_logT  ) / (n_T   - 1.0)
dye     = (max_ye     - min_ye    ) / (n_ye  - 1.0)

string = "mpirun -np 2 -env OMP_NUM_THREADS=16 ../../exe/nut_blackbody param.lua " + \
          str(min_logrho) + " " + str(max_logrho) + " " + str(rho0) + " " + str(n_rho) + " " + \
          str(min_logT  ) + " " + str(max_logT  ) + " " + str(T0  ) + " " + str(n_T  ) + " " + \
          str(min_ye    ) + " " + str(max_ye    ) + " " + str(ye0 ) + " " + str(n_ye ) + " " + \
          eosfile
print string
os.system(string)

results = np.loadtxt("results.dat",usecols=(0,2,3,4))
os.system("rm -f predicted.dat")
for i in range(len(results)):
    munue = results[i][0]
    rho   = results[i][1]
    T     = results[i][2]
    ye    = results[i][3]
    write_predicted(munue,rho,T,ye)

plot(rho0,T0,ye0)
