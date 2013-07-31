import numpy as np
import constants as pc
import pylab as py

T = 10.0e3
data = np.loadtxt('optical_I1.spec')
x = data[:,0]
y = data[:,1]
py.plot(x,y/np.average(y))

n = 2*pc.h*x**3/pc.c/(np.exp(pc.h*x/pc.k/T) - 1)
py.plot(x,n/np.average(n))
T = 5.0e3
n = 2*pc.h*x**3/pc.c/(np.exp(pc.h*x/pc.k/T) - 1)
py.plot(x,n/np.average(n))

py.show()
