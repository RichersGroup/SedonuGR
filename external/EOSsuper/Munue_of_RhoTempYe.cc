#include <cstdio>
#include <cstdlib>
#include <sstream>
#include <string>
#include "nuc_eos.hh"
#include <iostream>
#include <fstream>
#include <cmath>

using namespace nuc_eos;
const double k_MeV = 8.6173324e-11;   // boltzmann constant (Mev/K)

int main(int argc, char* argv[]) {

	if(argc!=5){
		std::cout << "Usage: Munue_of_RhoTempYe(MeV,incl. rest masses) EOSfile rho(g/ccm) T(MeV) Ye" << std::endl;
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
	std::cout << mue-muhat << std::endl;
	return 0;
}
