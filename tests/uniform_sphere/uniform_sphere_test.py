import numexpr as ne
import matplotlib as mpl
import scipy.integrate as integrate
mpl.use('Agg')
makeplot = False
if makeplot:
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
    mpl.rcParams['font.size'] = 16
    mpl.rcParams['font.family'] = 'serif'
    mpl.rcParams['xtick.major.size'] = 7
    mpl.rcParams['xtick.major.width'] = 2
    mpl.rcParams['xtick.major.pad'] = 8
    mpl.rcParams['xtick.minor.size'] = 4
    mpl.rcParams['xtick.minor.width'] = 2
    mpl.rcParams['ytick.major.size'] = 7
    mpl.rcParams['ytick.major.width'] = 2
    mpl.rcParams['ytick.minor.size'] = 4
    mpl.rcParams['ytick.minor.width'] = 2
    mpl.rcParams['axes.linewidth'] = 2
import h5py
import numpy as np
import sys
sys.path.insert(0, '../')
import tools

def chi(fluxfac,interpolator,e=None):

    if interpolator=="thick":
        chi = 1. / 3.
    if interpolator=="thin":
        chi = 1.
    if interpolator=="Kershaw":
        chi = ne.evaluate("(1. + 2.*fluxfac**2) / 3.")
    if interpolator=="Wilson":
        chi = ne.evaluate("1./3. - 1./3.*fluxfac + fluxfac**2")
    if interpolator=="Levermore":
        chi = ne.evaluate("(3 + 4.*fluxfac**2) / (5. + 2.*sqrt(4. - 3.*fluxfac**2))")
    if interpolator=="MEFD":
        fmax = ne.evaluate("1.-e")
        pmax = ne.evaluate("(2.*fmax*(1.-2.*e) + 1.) / 3.")
        pdiff = 1./3.
        x = ne.evaluate("fluxfac / fmax")
        x[np.where(x<0)] = 0
        x[np.where(x>1)] = 1.
        zeta = ne.evaluate("x**2*(3. - x + 3.*x**2)/5.")
        chi = ne.evaluate("(pmax-pdiff)*zeta + pdiff") # transforming from MEFD definition of chi to regular one
    if interpolator=="MEFDmp":
        chi = ne.evaluate("1./3. * (1.0 -2.*fluxfac + 4.*fluxfac**2)")
    if interpolator=="MEFDc":
        chi = ne.evaluate("1./3. + 2.*fluxfac**2/15. * (3. - fluxfac + 3.*fluxfac**2)")
    if interpolator=="ME":
        chi = ne.evaluate("1./3. + (2.*fluxfac**2)/15. * (3. - fluxfac + 3.*fluxfac**2)")
    if interpolator=="Janka1":
        a=0.5
        b=1.3064
        n=4.1342
        chi = ne.evaluate("1./3. * (1. + a*fluxfac**b + (2.-a)*fluxfac**n)")
    if interpolator=="Janka2":
        a=1.0
        b=1.345
        n=5.1717
        chi = ne.evaluate("1./3. * (1. + a*fluxfac**b + (2.-a)*fluxfac**n)")
    if interpolator=="mine":
        interpolator = ne.evaluate("fluxfac**5")
        chi = ne.evaluate("1./3.*(1.-interpolator) + 1.*interpolator")
    return chi




def chi_L(fluxfac, interpolator, e=None):
    if interpolator=="MEFDc":
        result = ne.evaluate("(2.*fluxfac**5 + 1.)/3.")
    elif interpolator=="MEFDmp":
        result = ne.evaluate("(1./3.)*(3. - 10.*fluxfac + 10.*fluxfac**2)")
    elif interpolator=="MEFD":
        fmax = ne.evaluate("1.-e")
        x = ne.evaluate("fluxfac/fmax")
        x[np.where(x<0)] = 0
        x[np.where(x>1)] = 1.
        zeta_MEFD = ne.evaluate("x**6")
        lmax = ne.evaluate("fmax * (1.-2.*e+2.*e**2)")
        ldiff = ne.evaluate("0.6*fluxfac")
        ldiffmax = ne.evaluate("0.6*fmax")
        lfree = fluxfac
        result = ne.evaluate("(2./3.)*zeta_MEFD*(lmax-ldiffmax)/(lfree-ldiff) + 1./3.")
        result[np.where(x==0)]=0
    else:
        result = chi(fluxfac,interpolator,e)
    
    return result

def p(f, interpolator):
    pfree = 1.
    pdiff = 1./3.
    X = chi(f, interpolator)
    return (3.*X-1)/2 * pfree + 3.*(1.-X)/2. * pdiff

def l(f, interpolator):
    lfree = f
    ldiff = 3.*f/5.
    X = chi_L(f, interpolator)
    return (3.*X-1)/2 * lfree + 3.*(1.-X)/2. * ldiff

# theoretical results
R=1
kappa=4
def mumin(r):
    if r>R:
        return np.sqrt(1-(R/r)**2)
    else:
        return -1
def g(r,mu):
    return np.sqrt(1-(r/R)**2 * (1-mu**2))
def s(r,mu):
    if(r<=R):
        return r*mu/R + g(r,mu)
    else:
        return 2.*g(r,mu)
def fdist(r,mu):
    return 1 - np.exp(-kappa*R*s(r,mu))

def e_theory(r):
    return integrate.quad(lambda mu: fdist(r,mu), mumin(r), 1)[0]
def f_theory(r):
    return integrate.quad(lambda mu: mu*fdist(r,mu), mumin(r), 1)[0]
def p_theory(r):
    return integrate.quad(lambda mu: mu**2*fdist(r,mu), mumin(r), 1)[0]
def l_theory(r):
    return integrate.quad(lambda mu: mu**3*fdist(r,mu), mumin(r), 1)[0]
    

avg_tolerance = 0.01
max_tolerance = 0.1

f = h5py.File("fluid_00001.h5","r")
r = np.array(f["axes/x0(cm)[mid]"])
E = np.array(f["distribution0(erg|ccm,tet)"])[:,0, 0]
F = np.array(f["distribution0(erg|ccm,tet)"])[:,0, 3]/E
P = np.array(f["distribution0(erg|ccm,tet)"])[:,0, 9]/E
L = np.array(f["distribution0(erg|ccm,tet)"])[:,0,19]/E

if makeplot:
    fig = plt.figure(figsize=(8,12))
    gs = gridspec.GridSpec(3, 1)
    axes = [plt.subplot(gsi) for gsi in gs]

isError = False

ClosureList = ["ME","MEFDmp","MEFDc","Levermore","Wilson","Kershaw","Janka1","Janka2"]

et = np.array([e_theory(rr) for rr in r])
ft = np.array([f_theory(rr) for rr in r]) / et
pt = np.array([p_theory(rr) for rr in r]) / et
lt = np.array([l_theory(rr) for rr in r]) / et

if makeplot:

    sedonu_linewidth = 3.0
    
    axes[0].set_ylabel(r"$f$")
    axes[0].plot(r, F, label="Sedonu",linewidth=sedonu_linewidth)
    axes[0].plot(r, ft, label="analytic",linewidth=3, color="k",linestyle="--")
    axes[0].legend(frameon=False)

    axes[1].set_ylabel(r"$p$")
    axes[1].plot(r, P, linewidth=sedonu_linewidth)
    for pClosure in ClosureList:
        axes[1].plot(r, p(F,pClosure), label=pClosure, alpha=1.0)
    axes[1].plot(r, pt,linewidth=3, color="k",linestyle="--")
    axes[1].legend(frameon=False, ncol=2, loc=(.35,1.15), handlelength=1.5)
    
    axes[2].set_ylabel(r"$l$")
    axes[2].plot(r, L, linewidth=sedonu_linewidth)
    for lClosure in ClosureList:
        axes[2].plot(r, l(F,lClosure), alpha=1)
    axes[2].plot(r, lt,linewidth=3, color="k",linestyle="--")

    for ax in axes:
        ax.set_xlim(0.5,2)
    axes[0].set_ylim(0,1)
    axes[1].set_ylim(.2,.95)
    axes[2].set_ylim(0,.95)
    
def check_error(error,name):
    avgerror = np.abs(np.mean(error))
    maxerror = np.max(np.abs(error))
    print(name+" density avgerror=",avgerror," maxerror=",maxerror)
    if avgerror>avg_tolerance or maxerror>max_tolerance:
        isError = True

check_error(F-ft,"Fluxfac")
check_error(P-pt,"Pressure")
check_error(L-lt,"Rank3")
    
if makeplot:
    for ax in axes:
        ax.minorticks_on()
        ax.tick_params(axis='both',which='both',direction='in',top=True,right=True)
    
    axes[2].set_xlabel("$r$ (cm)")
    axes[0].xaxis.set_ticklabels([])
    axes[1].xaxis.set_ticklabels([])

    plt.subplots_adjust(wspace=0, hspace=0)
    plt.savefig("uniform_sphere.pdf", bbox_inches="tight")

exit()

if isError:
    raise Exception("temperature_redshift results are outside of the tolerance.")
else:
    print("SUCCESS")
