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
tolerance = 6e-2
passing = True

def calc_error(string, data, theory, passing):
    chi2 = np.sqrt((data-theory)**2)
    error = chi2.sum()/np.size(data)
    pass_this_test = (abs(error)<tolerance)
    print(string, error, pass_this_test)
#    if(not pass_this_test):
#        print(theory, data)
    return passing and pass_this_test

f = h5py.File("fluid_00001.h5","r")
edens = np.array(f["distribution0(erg|ccm,tet)"])
shape = np.shape(edens)
ncells = shape[0]*shape[1]*shape[2]
#edens_avg = edens.sum(axis=(0,1,2)) / ncells
print("=== ERRORS ===")
passing = calc_error("Fx/E:" ,edens[:,:,:,:,1]/edens[:,:,:,:,0], 0    ,passing)
passing = calc_error("Fy/E:" ,edens[:,:,:,:,2]/edens[:,:,:,:,0], 0    ,passing)
passing = calc_error("Fz/E:" ,edens[:,:,:,:,3]/edens[:,:,:,:,0], 0    ,passing)
passing = calc_error("Pxx/E:",edens[:,:,:,:,4]/edens[:,:,:,:,0], 1./3.,passing)
passing = calc_error("Pxy/E:",edens[:,:,:,:,5]/edens[:,:,:,:,0], 0    ,passing)
passing = calc_error("Pxz/E:",edens[:,:,:,:,6]/edens[:,:,:,:,0], 0    ,passing)
passing = calc_error("Pyy/E:",edens[:,:,:,:,7]/edens[:,:,:,:,0], 1./3.,passing)
passing = calc_error("Pyz/E:",edens[:,:,:,:,8]/edens[:,:,:,:,0], 0    ,passing)
passing = calc_error("Pzz/E:",edens[:,:,:,:,9]/edens[:,:,:,:,0], 1./3.,passing)
#edens_theory = np.array([tools.edens(T_gas[i]*tools.MeV, munue*tools.MeV) for i in range(len(T_gas))])

force_abs = np.array(f["four-force[abs](erg|ccm|s,tet)"])
force_emit = np.array(f["four-force[emit](erg|ccm|s,tet)"])
passing = calc_error("force_x:",(force_abs[:,:,:,0]+force_emit[:,:,:,0])/force_emit[:,:,:,0], 0,passing)
passing = calc_error("force_y:",(force_abs[:,:,:,1]+force_emit[:,:,:,1])/force_emit[:,:,:,1], 0,passing)
passing = calc_error("force_z:",(force_abs[:,:,:,2]+force_emit[:,:,:,2])/force_emit[:,:,:,2], 0,passing)
passing = calc_error("force_t:",(force_abs[:,:,:,3]+force_emit[:,:,:,3])/force_emit[:,:,:,3], 0,passing)

edens_theory = tools.edens(T, munue)
passing = calc_error("edens:",edens[:,:,:,:,0].sum(axis=3)/edens_theory, 1,passing)

if(passing):
    print("SUCCESS")
else:
    raise Exception("3Dbox results are outside of the tolerance.")
