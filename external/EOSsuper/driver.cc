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

double beta_eq_ye(double rho, double T){
  double solve_for_munu = 0;
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
  //std::cout << iter << std::endl;
  //if(iter==maxiter) std::cout << "WARNING - could not converge on a ye" << std::endl;
  return xye;
}



int main(int argc, char* argv[]) {

  const double k_MeV = 8.6173324e-11;   // boltzmann constant (Mev/K)
  nuc_eos_C_ReadTable("../../external/tables/LS220.h5");
  double x, rho, xrho, xtemp, trash;
  
  // set up the files
  std::fstream fin(argv[1]);
  std::string line;
  std::stringstream ss;
  getline(fin,line); //ignore the first line
  while(getline(fin,line)){
    ss.str(line);
    ss >> x;
    ss >> rho;
    xrho = rho*RHOGF;
    ss >> xtemp;
    xtemp *= k_MeV;
    
    std::cout << x << "\t" << rho << "\t" << beta_eq_ye(xrho,xtemp) << std::endl;
  }

  return 0;
}
