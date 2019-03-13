import h5py
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import sys
sys.path.insert(0, '../')
import tools
import matplotlib.gridspec as gridspec

munue = 10. * tools.MeV
T = 10. * tools.MeV
tolerance = 5e-3
passing = True

def calc_error(string, data, theory, passing):
    error = (data-theory).sum()/np.size(data)
    pass_this_test = (abs(error)<tolerance)
    print(string, error)
    return passing and pass_this_test

f = h5py.File("fluid_00001.h5","r")
edens = np.array(f["distribution0(erg|ccm,tet)"])
shape = np.shape(edens)
ncells = shape[0]*shape[1]*shape[2]
edens_avg = edens.sum(axis=(0,1,2)) / ncells
print("=== ERRORS ===")
passing = calc_error("Fx/E:",edens_avg[:,1]/edens_avg[:,0], 0    ,passing)
passing = calc_error("Fy/E:",edens_avg[:,2]/edens_avg[:,0], 0    ,passing)
passing = calc_error("Fz/E:",edens_avg[:,3]/edens_avg[:,0], 0    ,passing)
passing = calc_error("Pxx/E:",edens_avg[:,4]/edens_avg[:,0], 1./3.,passing)
passing = calc_error("Pxy/E:",edens_avg[:,5]/edens_avg[:,0], 0    ,passing)
passing = calc_error("Pxz/E:",edens_avg[:,6]/edens_avg[:,0], 0    ,passing)
passing = calc_error("Pyy/E:",edens_avg[:,7]/edens_avg[:,0], 1./3.,passing)
passing = calc_error("Pyz/E:",edens_avg[:,8]/edens_avg[:,0], 0    ,passing)
passing = calc_error("Pzz/E:",edens_avg[:,9]/edens_avg[:,0], 1./3.,passing)
#edens_theory = np.array([tools.edens(T_gas[i]*tools.MeV, munue*tools.MeV) for i in range(len(T_gas))])

force_abs = np.array(f["four-force[abs](erg|ccm|s,tet)"]).sum(axis=(0,1,2))/ncells
force_emit = np.array(f["four-force[emit](erg|ccm|s,tet)"]).sum(axis=(0,1,2))/ncells
passing = calc_error("force_x:",(force_abs[0]+force_emit[0])/force_emit[0], 0,passing)
passing = calc_error("force_y:",(force_abs[1]+force_emit[1])/force_emit[1], 0,passing)
passing = calc_error("force_z:",(force_abs[2]+force_emit[2])/force_emit[2], 0,passing)
passing = calc_error("force_t:",(force_abs[3]+force_emit[3])/force_emit[3], 0,passing)

edens_theory = tools.edens(T, munue)
passing = calc_error("edens:",edens_avg[:,0].sum()/edens_theory, 1,passing)

if(passing):
    print("SUCCESS")
else:
    raise Exception("3Dbox results are outside of the tolerance.")
