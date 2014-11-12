#!/usr/bin/python
import sys
from numpy import *
import numpy as np

if len(sys.argv)!=2:
    print("Usage: python convert_outer_to_center.py input_file")
    sys.exit()

data = loadtxt(sys.argv[1])
radii = data[:,0]

centers = np.zeros(len(radii))
for i in range(0,len(radii)):
    if(i==0):
        print  0
    else:
        print 0.5*(radii[i]+radii[i-1])
