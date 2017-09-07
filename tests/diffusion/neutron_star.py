import sys

dx0 = 0.005
Rmin = 1
Rmax = 1.1

Rlist = []
R=Rmin
while R<Rmax:
    R += dx0 * R
    Rlist.append(R)

print('1D_sphere',len(Rlist),Rmin)
for i  in range(len(Rlist)):
    print(Rlist[i], 1, 0, 0, 0, 0, 0)
