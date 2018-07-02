import os
from math import *
from scipy.integrate import quad
import numpy as np

# INPUTS
tolerance = 0.05
eosfile   = "SFHo.h5" #ignored if compiled for helmholtz eos

min_logrho = 6  #g/ccm
max_logrho = 15 #g/ccm
n_rho = 5
center_logrho = 10

min_logT = 0 #MeV
max_logT = 1.5  #MeV
n_T = 5
center_logT = 0.5

min_ye = 0.05
max_ye = 0.55
n_ye = 5
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

def predicted_edens(T,mu):
    return edens(T,mu) + edens(T,-mu) + 4*edens(T,0)
def calculated_edens(results):
    edens = []
    for i in range(len(results)):
        edens.append(results[i][4]+results[i][5]+results[i][6])
    return edens

def write_predicted(mu,rho,T,Ye):
    etot = predicted_edens(T,mu)
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

string = "OMP_NUM_THREADS=16 mpirun -np 1 ../../exe/blackbody_test param.lua " + \
          str(min_logrho) + " " + str(max_logrho) + " " + str(rho0) + " " + str(n_rho) + " " + \
          str(min_logT  ) + " " + str(max_logT  ) + " " + str(T0  ) + " " + str(n_T  ) + " " + \
          str(min_ye    ) + " " + str(max_ye    ) + " " + str(ye0 ) + " " + str(n_ye ) + " " + \
          eosfile
print(string)
os.system(string)

results = np.loadtxt("results.dat",usecols=(0,1,2,3,4,5,6))
os.system("rm -f predicted.dat")
for i in range(len(results)):
    munue = results[i][3]
    rho   = results[i][0]
    T     = results[i][1]
    ye    = results[i][2]
    write_predicted(munue,rho,T,ye)
    print(rho,T,ye)
    print(predicted_edens(T,munue), calculated_edens(results)[i])
    error = abs(predicted_edens(T,munue) - calculated_edens(results)[i]) / (predicted_edens(T,munue) + calculated_edens(results)[i])
    if error>tolerance:
        print(predicted_edens(T,munue), calculated_edens(results)[i], error)
        raise Exception("blackbody results are outside of the tolerance.")
    else:
        print("SUCCESS: error=",error)


    
#plot(rho0,T0,ye0)
