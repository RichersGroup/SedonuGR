import numpy as np
import constants as pc
import pylab as py

data = np.loadtxt('ray_00001')
r = data[:,0]
y = data[:,2]

L  = 1e43
r0 = 2e15

w = 0.5*(1 - (1 - r0**2/r**2)**0.5)

T = (L*w/(4.0*pc.pi*r0**2)/pc.sb)**0.25

#T = 7600*(r/r0)**(-0.5)

py.plot(r,y)
py.plot(r,T)
#py.yscale('log')
py.show()
