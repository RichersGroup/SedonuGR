#include <cstdio>
#include <cstdlib>
#include <sstream>
#include <string>
#include "nuc_eos.hh"
#include <iostream>
#include <fstream>
#include <cmath>

using namespace nuc_eos;
const double me = 0.510998910; // in MeV
const double mp = 938.272046;  // in MeV
const double mn = 939.565378;  // in MeV
const double k_MeV = 8.6173324e-11;   // boltzmann constant (Mev/K)


int main(int argc, char* argv[]) {

	if(argc!=5){
		std::cout << "Usage: Ye_of_RhoTempMunue EOSfile rho(g/ccm) T(MeV) Munue(MeV,incl. rest masses)" << std::endl;
		exit(1);
	}

	//read in parameters
	nuc_eos_C_ReadTable(argv[1]);
	double rho = atof(argv[2])*RHOGF; // input g/ccm --> code units
	double T = atof(argv[3]);         // MeV
	double mu = atof(argv[4]);        // MeV
	//std::cout << "rho: " << rho/RHOGF << std::endl;
	//std::cout << "T: "<< T << std::endl;
	//std::cout << "Ye: " << Ye << std::endl;

	// set up the solver
	double solve_for_munu = mu + (mn-mp); // eos does not include baryon masses
	double ye_min = 0.05;
	double ye_max = 0.55;
	double xye = 0.5;
	double xeps, xprs, xent, xcs2, xdedt, xdpderho, xdpdrhoe, xmunu=10;
	int keyerr, anyerr;
	int n=1;
	int iter=0;
	int maxiter = 1000;
	while(fabs(xmunu-solve_for_munu) > 1e-10 && iter<maxiter){
		xye = 0.5*(ye_min + ye_max);
		nuc_eos_m_kt1_short(&n,&rho,&T,&xye,
				&xeps,&xprs,&xent,&xcs2,&xdedt,
				&xdpderho,&xdpdrhoe,&xmunu,
				&keyerr,&anyerr);
		if(xmunu-solve_for_munu<0) ye_min = xye;
		else ye_max = xye;
		iter++;
	}
	std::cout << xye << std::endl;

	return 0;
}
