#include <cstdio>
#include <cstdlib>
#include <sstream>
#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include "../external/EOSsuper/nuc_eos.hh"
#include "global_options.h"

void EOSsuper_initialize(char* filename){
	nuc_eos::nuc_eos_C_ReadTable(filename);
}

double EOSsuper_munue(const double rho_cgs /*g/ccm*/, const double T_cgs/*K*/, const double Ye){
	double T = T_cgs * pc::k_MeV;
	double rho = rho_cgs * RHOGF;

	// various parameters and output
	int n=1;
	double xeps, xprs, xent, xcs2, xdedt, xdpderho, xdpdrhoe;
	double xa, xh, xn, xp, abar, zbar;
	double mue, mun, mup, muhat;
	int keyerr, anyerr;

	// call the equation of state
	nuc_eos::nuc_eos_m_kt1_full(&n,&rho,&T,&Ye,
			   &xeps,&xprs,&xent,&xcs2,&xdedt,
			   &xdpderho,&xdpdrhoe,&xa,&xh,&xn,&xp,&abar,&zbar,
			   &mue,&mun,&mup,&muhat,
			   &keyerr,&anyerr);

	// output the answer
	return (mue-muhat) * pc::MeV_to_ergs; // ergs
}
