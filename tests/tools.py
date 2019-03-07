from scipy.integrate import quad
import numpy as np
import math

# Constants
h = 6.6260755e-27 #erg s
c = 2.99792458e10   #cm/s
MeV = 1.60217646e-6  #erg/MeV
k_b = 1.380658e-16 # erg/K

# x = h nu / kT
# y = mu / kT
def FermiDirac(x,y):
    return 1. / ( math.exp(x-y) + 1. )
def numberBlackbody(x,y):
    return x**2 * FermiDirac(x,y)
def energyBlackbody(x,y):
    return x**3 * FermiDirac(x,y)

# kT and mu in units of ergs
def ndens_bin(kT,mu,x0,x1):
    y = mu/kT
    return 4.*np.pi * (kT/(c*h))**3 * \
        quad(numberBlackbody,x0,x1,args=(y))[0]
def edens_bin(kT,mu,x0,x1):
    y = mu/kT
    return 4.*np.pi * kT * (kT/(c*h))**3 * \
        quad(energyBlackbody,x0,x1,args=(y))[0]

def ndens(kT,mu):
    if(mu/kT > 10): max_x = 5*mu/kT
    else: max_x = 20
    return ndens_bin(kT,mu,0,max_x)
def edens(kT,mu):
    if(mu/kT > 10): max_x = 5*mu/kT
    else: max_x = 20
    return edens_bin(kT,mu,0,max_x)
