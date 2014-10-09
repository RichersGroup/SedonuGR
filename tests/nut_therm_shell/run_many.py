import os
import subprocess
from math import *
from scipy.integrate import quad
from numpy import Inf

# INPUTS
eosfile   = "/data/tables/EOS/HShen.h5"

min_logrho = 8  #g/ccm
max_logrho = 15 #g/ccm
n_rho = 100

min_logT = -1 #MeV
max_logT = 2  #MeV
n_T = 100

min_ye = 0.05
max_ye = 0.55
n_ye = 100

# Constants
h = 6.626075518e-27 #erg s
c = 2.99792458e10   #cm/s
MeV = 1.60217657e-6  #erg/MeV


def munue(rho,T,Ye,eosfile):
    command = ['../../external/EOSsuper/Munue_of_RhoTempYe',eosfile,str(rho),str(T),str(Ye)]
    p = subprocess.Popen(command, stdout=subprocess.PIPE)
    p.wait()
    return float(p.stdout.read().split('\n')[-2])

def write_results(rho,T,Ye,munue,filename):
    string = "'"+str(rho)+" "+str(T)+" "+str(Ye)+" "+str(munue)+" '"
    os.system("echo -n "+string+" >> "+filename)
    os.system("tail -n1 fluid_00005.dat >> "+filename)

def plot():
    string = "sed " + \
             "-e 's/RHO_HERE/"    + str(rho0)         + "/g' " + \
             "-e 's/TEMP_HERE/"   + str(T0)           + "/g' " + \
             "-e 's/YE_HERE/"     + str(ye0)          + "/g' " + \
             "template.gnuplot > compare.gnuplot"
    os.system(string)
    os.system("gnuplot compare.gnuplot")

def run_test(rho,T,Ye,solve_T,solve_ye,filename):
    print "Currently running: rho="+str(rho)+"g/ccm T="+str(T)+"MeV Ye="+str(Ye)
    os.system("rm -f fluid_*.dat")

    mu = munue(rho,T,Ye,eosfile)
    string = "sed " + \
             "-e 's/SOLVE_T_HERE/"    + str(solve_T)    + "/g' " + \
             "-e 's/SOLVE_YE_HERE/"   + str(solve_ye)   + "/g' " + \
             "-e 's/TEMP_HERE/"       + str(T)          + "/g' " + \
             "-e 's/MUNU_HERE/"       + str(mu)         + "/g' " + \
             "template.lua > param.lua"
    os.system(string)
    string = "sed " + \
             "-e 's/RHO_HERE/" + str(rho) + "/g' " + \
             "template.py > shell.py"
    os.system(string)
    os.system("python shell.py > shell.mod")
    os.system("export OMP_NUM_THREADS=5")
    os.system("mpirun -np 16 ./gomc")

    write_results(rho,T,Ye,mu,filename)




### BEGIN SCRIPT ###

file_T  = "results_T.dat"
file_ye = "results_ye.dat"
file_both = "results.dat"
os.system("rm -f "+file_T)
os.system("rm -f "+file_ye)
os.system("rm -f "+file_both)

rho0 = 10**(0.5*(min_logrho + max_logrho))
T0   = 10**(0.5*(min_logT   + max_logT  ))
ye0  =      0.5*(min_ye     + max_ye    )
dlogrho = (max_logrho - min_logrho) / (n_rho - 1.0)
dlogT   = (max_logT   - min_logT  ) / (n_T   - 1.0)
dye     = (max_ye     - min_ye    ) / (n_ye  - 1.0)

for i in range (4, n_rho):
    logrho = min_logrho + i*dlogrho
    run_test(10**logrho,T0,ye0,1,0,file_T)
    run_test(10**logrho,T0,ye0,0,1,file_ye)
    run_test(10**logrho,T0,ye0,1,1,file_both)

for i in range (0, n_T):
    logT = min_logT + i*dlogT
    run_test(rho0,10**logT,ye0,1,0,file_T)
    run_test(rho0,10**logT,ye0,0,1,file_ye)
    run_test(rho0,10**logT,ye0,1,1,file_both)

for i in range (0, n_ye):
    ye = min_ye + i*dye
    run_test(rho0,T0,ye,1,0,file_T)
    run_test(rho0,T0,ye,0,1,file_ye)
    run_test(rho0,T0,ye,1,1,file_both)

plot()
