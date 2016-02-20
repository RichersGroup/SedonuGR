#include <cstdio>
#include <cstdlib>
#include <sstream>
#include <string>
#include "nuc_eos.hh"
#include <iostream>
#include <fstream>
#include <cmath>

using namespace nuc_eos;
using namespace std;
const double k_MeV = 8.6173324e-11;   // boltzmann constant (Mev/K)

int main(int argc, char* argv[]) {

	if(argc!=5){
		std::cout << "Usage: Mue_of_RhoTempYe(MeV,incl. rest masses) EOSfile rho(g/ccm) T(MeV) Ye" << std::endl;
		exit(1);
	}

	//read in parameters
	nuc_eos_C_ReadTable(argv[1]);
	double rho = atof(argv[2])*RHOGF;
	double T = atof(argv[3]);
	double Ye = atof(argv[4]);

	// various parameters and output
	int n=1;
	double xeps, xprs, xent, xcs2, xdedt, xdpderho, xdpdrhoe;
	double xa, xh, xn, xp, abar, zbar;
	double mue, mun, mup, muhat;
	int keyerr, anyerr;

	// call the equation of state
	nuc_eos_m_kt1_full(&n,&rho,&T,&Ye,
			   &xeps,&xprs,&xent,&xcs2,&xdedt,
			   &xdpderho,&xdpdrhoe,&xa,&xh,&xn,&xp,&abar,&zbar,
			   &mue,&mun,&mup,&muhat,
			   &keyerr,&anyerr);

	// output the answer
	cout << "eps: " << xeps << " erg/g" << endl;
	cout << "pressure: " << xprs << " dyn/cm^2" << endl;
	cout << "entropy: " << xent << " kB/baryon" << endl;
	cout << "sound speed^2: " << xcs2 << " (cm/s)^2" << endl;
	cout << "de/dT: " << xdedt << "" << endl;
	cout << "dP/derho: " << xdpderho << "" << endl;
	cout << "dP/drhoe: " << xdpdrhoe << "" << endl;
	cout << "X_alpha: " << xa << endl;
	cout << "X_hydrogen: " << xh << endl;
	cout << "X_neutron: " << xn << endl;
	cout << "X_proton: " << xp << endl;
	cout << "Abar: " << abar << endl;
	cout << "Zbar: " << zbar << endl;
	cout << "mu_e: " << mue << " MeV" << endl;
	cout << "mu_n: " << mun << " MeV" << endl;
	cout << "mu_p: " << mup << " MeV" << endl;
	cout << "muhat: " << muhat << " MeV" << endl;
	cout << "keyerr: " << keyerr << endl;
	cout << "anyerr: " << anyerr << endl;
	cout << "munue: " << mue-muhat << endl;
	return 0;
}
