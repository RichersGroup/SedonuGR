import numpy as np
import math
import os
import sys
import glob
import shutil
import h5py

h = 6.6260755e-27 # erg s
MeV = 1.60218e-6 # erg
c = 2.99792458e10 # cm/s
kB=1.380658e-16 # erg/K

ix=103
iy=103
iz=0

# use file 0 as template
infilename = sys.argv[1]
infile = h5py.File(infilename,"r")

# print point properties
print("ix=",ix,"iy=",iy,"iz=",iz)
print("rho=",infile["rho(g|ccm,tet)"][ix,iy,iz]/1e12,"*1e12")
print("Ye=",infile["Ye"][ix,iy,iz])
print("T=",infile["T_gas(K,tet)"][ix,iy,iz]*kB/MeV," MeV")
print("sqrtdetg=",infile["sqrtdetg3"][ix,iy,iz])
print("v=",infile["threevelocity(cm|s)"][ix,iy,iz]/c)

# get the energy grid and phase space volumes (erg)
#Emid  = np.array(infile["axes/frequency(Hz)[mid]" ]) * h # erg
nuedge = np.array(infile["axes/frequency(Hz)[edge]"]) # 1/s
dnu4 = np.array([nuedge[i+1]**4 - nuedge[i]**4 for i in range(len(nuedge)-1)])
#phasevol = Emid * np.array([nuedge[i+1]**3-nuedge[i]**3 for i in range(len(Emid))]) * 4.*np.pi
maxEdens = np.pi*h/c**3 * dnu4

# prepare datasets to combine
distribution0 = np.array(infile["distribution0(erg|ccm,tet)"])[ix,iy,iz,:,0:4] / maxEdens[:,np.newaxis]
distribution1 = np.array(infile["distribution1(erg|ccm,tet)"])[ix,iy,iz,:,0:4] / maxEdens[:,np.newaxis]
distribution2 = np.array(infile["distribution2(erg|ccm,tet)"])[ix,iy,iz,:,0:4] / maxEdens[:,np.newaxis]/4.

J = np.array([distribution0[:,0  ],distribution1[:,0  ],distribution2[:,0  ]])
H = np.array([distribution0[:,1:4],distribution1[:,1:4],distribution2[:,1:4]])

for s in range(3):
    print(s, J[s])
    #print(s, np.min(H[s]),np.max(H[s]))
    
print(np.shape(distribution0))
print(np.where(J[0]<0))
