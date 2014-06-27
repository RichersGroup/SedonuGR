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
		std::cout << "Usage: Mue_of_RhoTempYe EOSfile rho(g/ccm) T(MeV) Ye" << std::endl;
		exit(1);
	}

	//read in parameters
	nuc_eos_C_ReadTable(argv[1]);
	double rho = atof(argv[2])*RHOGF;
	double T = atof(argv[3]);
	double Ye = atof(argv[4]);
	std::cout << "rho: " << rho/RHOGF << std::endl;
	std::cout << "T: "<< T << std::endl;
	std::cout << "Ye: " << Ye << std::endl;

	// various parameters and output
	int n=1;
	double xeps, xprs, xent, xcs2, xdedt, xdpderho, xdpdrhoe, xmunu;
	int keyerr, anyerr;

	// call the equation of state
	nuc_eos_m_kt1_short(&n,&rho,&T,&Ye,
			&xeps,&xprs,&xent,&xcs2,&xdedt,
			&xdpderho,&xdpdrhoe,&xmunu,
			&keyerr,&anyerr);

	// output the answer
	std::cout << xmunu << std::endl;
	return 0;
}
