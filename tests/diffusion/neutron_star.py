import sys

dx0 = 0.005
Rmin = 100
Rmax = 100.1
alpha = 1.0
X = 1.0

Rlist = []
R=Rmin
while R<Rmax:
    R += dx0
    Rlist.append(R)

print('1D_sphere',len(Rlist),Rmin)
for i  in range(len(Rlist)):
    print(Rlist[i], 1, 0, 0, 0, alpha, X)
