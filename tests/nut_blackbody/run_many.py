import os
import subprocess
from math import *
from scipy.integrate import quad
from numpy import Inf

# INPUTS
eosfile   = "/data/tables/EOS/HShen.h5"
nulibfile = "\"\/data\/tables\/NuLib\/NuLib_Helmholtz_noscat.h5\"" # requires backslashes b/c it gets put within regex

min_logrho = 8  #g/ccm
max_logrho = 15 #g/ccm
n_rho = 50
center_logrho = 10

min_logT = -1 #MeV
max_logT = 2  #MeV
n_T = 50
center_logT = 0.5;

min_ye = 0.05
max_ye = 0.55
n_ye = 50
center_ye = 0.3

nparticles = 1e6


# Constants
h = 6.626075518e-27 #erg s
c = 2.99792458e10   #cm/s
MeV = 1.60217657e-6  #erg/MeV


def munue(rho,T,Ye,eosfile):
    command = ['../../external/EOSsuper/Munue_of_RhoTempYe',eosfile,str(rho),str(T),str(Ye)]
    p = subprocess.Popen(command, stdout=subprocess.PIPE)
    p.wait()
    return float(p.stdout.read().split('\n')[-2])

def BB(x,T,mu):
    # x = h nu / kT
    # T and mu are in same units
    return x**3 / ( exp(x - mu/T) + 1 )

def edens(T,mu):
    # mu and T are in the same units
    if(mu/T > 10): max_x = 3*mu/T
    else: max_x = 20
    return 4*pi*h*c * (T*MeV/c/h)**4 * quad(BB,0.0,max_x,args=(T,mu))[0]

def write_results():
    os.system("tail -n1 fluid_00001.dat >> results.dat")

def write_predicted(rho,T,Ye):
    mu = munue(rho,T,Ye,eosfile)
    etot = edens(T,mu) + edens(T,-mu) + 4*edens(T,0)
    string = str(etot)+" "+str(rho)+" "+str(T)+" "+str(Ye)+" "+str(mu)
    os.system("echo "+string+" >> predicted.dat")

def plot(rho0,T0,ye0):
    string = "sed " + \
             "-e 's/RHO_HERE/"    + str(rho0)         + "/g' " + \
             "-e 's/TEMP_HERE/"   + str(T0)           + "/g' " + \
             "-e 's/YE_HERE/"     + str(ye0)          + "/g' " + \
             "template.gnuplot > compare.gnuplot"
    os.system(string)
    os.system("gnuplot compare.gnuplot")

def run_test(rho,T,Ye):
    print "Currently running: rho="+str(rho)+"g/ccm T="+str(T)+"MeV Ye="+str(Ye)

    string = "sed " + \
             "-e 's/RHO_HERE/"        + str(rho)        + "/g' " + \
             "-e 's/TEMP_HERE/"       + str(T)          + "/g' " + \
             "-e 's/YE_HERE/"         + str(Ye)         + "/g' " + \
             "-e 's/NULIB_HERE/"      + nulibfile       + "/g' " + \
             "-e 's/NPARTICLES_HERE/" + str(nparticles) + "/g' " + \
             "template.lua > param.lua"
    os.system(string)
    os.system("mpirun -np 16 -env OMP_NUM_THREADS 2 ./gomc")

    write_results()
    write_predicted(rho,T,Ye)




### BEGIN SCRIPT ###

os.system("rm results.dat")
os.system("rm predicted.dat")

rho0 = 10**(center_logrho)
T0   = 10**(center_logT  )
ye0  =      center_ye
dlogrho = (max_logrho - min_logrho) / (n_rho - 1.0)
dlogT   = (max_logT   - min_logT  ) / (n_T   - 1.0)
dye     = (max_ye     - min_ye    ) / (n_ye  - 1.0)

for i in range (0, n_rho):
    logrho = min_logrho + i*dlogrho
    run_test(10**logrho,T0,ye0)

for i in range (0, n_T):
    logT = min_logT + i*dlogT
    run_test(rho0,10**logT,ye0)

for i in range (0, n_ye):
    ye = min_ye + i*dye
    run_test(rho0,T0,ye)

plot(rho0,T0,ye0)
