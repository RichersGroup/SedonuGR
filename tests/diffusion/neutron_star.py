import sys

dx0 = 0.01
Rmax = 100
Rmin = 1

Rlist = []
R=Rmin
while R<Rmax:
    R += dx0 * R
    Rlist.append(R)

print '1D_sphere',len(Rlist),Rmin
for i  in range(len(Rlist)):
    print Rlist[i], 0, 0, 0
