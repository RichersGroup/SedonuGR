import numpy as np
import sys
import subprocess
import os
from multiprocessing import Pool

nfiles = 3

qsetlist = [
    ['Pxx','Pyy','Pzz'],
    ['Pxy','Pxz','Pyz'],
    ['Lxxx','Lyyy','Lzzz'],
    ['Lxxy','Lxxz','Lxyy','Lxzz','Lxyz','Lyyz','Lyzz']
    ]
interpolatorlist = [
    'ME',
    'Kershaw',
    'Wilson',
    'Levermore',
    'MEFD_maxpack'
    ]
for chi in np.linspace(1./3., 1.0, num=20):
    interpolatorlist.append(str(chi))

def work_for_one_interpolator(interpolator):
    for qset in qsetlist:
        print(interpolator+" "+str(qset))
        commandlist = ["python"]
        commandlist.append("/home/srichers/software/sedonu/useful_scripts/analysis/closure_compare.py")
        commandlist.append(interpolator)
        for q in qset:
            commandlist.append(q)
            
        try:
            subprocess.call(commandlist)
        except KeyboardInterrupt:
            exit(0)

with Pool(7) as pool:
    pool.map(work_for_one_interpolator, interpolatorlist)
