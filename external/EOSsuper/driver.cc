#include <cstdio>
#include <cstdlib>
#include "nuc_eos.hh"
#include <iostream>
#include <fstream>
#include <cmath>

using namespace nuc_eos;
const double me = 0.510998910; // in MeV
const double mp = 938.272046;  // in MeV
const double mn = 939.565378;  // in MeV

double beta_eq_ye(double rho, double T){
  double solve_for_munu = -(mp+me-mn);
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
  std::cout << iter << std::endl;
  if(iter==maxiter) std::cout << "WARNING - could not converge on a ye" << std::endl;
  return xye;
}



int main(int argc, char* argv[]) {


  nuc_eos_C_ReadTable("/home/sherwood/software/sedona_git/external/tables/LS220.h5");
  //printf("energy_shift: %15.6E\n",energy_shift);
  //printf("nrho: %d\n",nrho);
  //printf("ntemp: %d\n",ntemp);
  //printf("nye: %d\n",nye);

  double xrho = 3.0e14 * RHOGF;
  double xtemp = 0.8;
  double xye = 0.5;
  double xeps = 0.0;
  double xprs = 0.0;
  int keyerr = 0;
  int anyerr = 0;
  int n = 1;
  
  nuc_eos_m_kt1_press_eps(&n,&xrho,&xtemp,&xye,&xeps,&xprs,&keyerr,&anyerr);

  fprintf(stderr,"*********************\n");
  fprintf(stderr,"kt1_press:\n");
  fprintf(stderr,"rho, T, ye: %15.6E %15.6E %15.6E\n",xrho,xtemp,xye);
  fprintf(stderr,"press, eps, keyerr, anyerr: %15.6E %15.6E %d %d\n",xprs,xeps,keyerr,anyerr);

  // now see what we get
  const double prec = 1.0e-12;
  //  xeps = 4.0e-4;
  nuc_eos_m_kt0_press(&n,&xrho,&xtemp,&xye,&xeps,&xprs,
		      &prec,&keyerr,&anyerr);

  fprintf(stderr,"*********************\n");
  fprintf(stderr,"kt0_press:\n");
  fprintf(stderr,"rho, T, ye: %15.6E %15.6E %15.6E\n",xrho,xtemp,xye);
  fprintf(stderr,"press, eps, keyerr, anyerr: %15.6E %15.6E %d %d\n",xprs,xeps,keyerr,anyerr);

  fprintf(stderr,"*********************\n");
  fprintf(stderr,"Short EOS Call:\n");
  xrho = 1.0e11 * RHOGF;
  xtemp = 2.0;
  xye = 0.1028415;

  // we assume we know our temperature
  // declare additional vars
  double xent,xcs2,xdedt,xdpderho,xdpdrhoe,xmunu;
  nuc_eos_m_kt1_short(&n,&xrho,&xtemp,&xye,
		      &xeps,&xprs,&xent,&xcs2,&xdedt,
		      &xdpderho,&xdpdrhoe,&xmunu,
		      &keyerr,&anyerr);
  fprintf(stderr,"rho, T, ye: %15.6E %15.6E %15.6E\n",xrho,xtemp,xye);
  fprintf(stderr,"press, eps, keyerr, anyerr: %15.6E %15.6E %d %d\n",xprs,xeps,keyerr,anyerr);
  fprintf(stderr,"xent: %15.6E xcs2: %15.6E xdedt: %15.6E\n",xent,xcs2,xdedt);
  fprintf(stderr,"xdpderho: %15.6E xdpdrhoe: %15.6E \n",xdpderho,xdpdrhoe);
  fprintf(stderr,"xmunu: %15.6E\n",xmunu);
  

  double xxa,xxh,xxn,xxp;
  double xabar,xzbar,xmue,xmun,xmup,xmuhat;
  nuc_eos_m_kt1_full(&n,&xrho,&xtemp,&xye,&xeps,
		     &xprs,&xent,&xcs2,&xdedt,&xdpderho,&xdpdrhoe,
		     &xxa,&xxh,&xxn,&xxp,&xabar,&xzbar,
		     &xmue,&xmun,&xmup,&xmuhat,&keyerr,&anyerr);

  fprintf(stderr,"*********************\n");
  fprintf(stderr,"%15.6E %15.6E %15.6E\n",xrho,xtemp,xye);
  fprintf(stderr,"%15.6E %15.6E %d %d\n",xprs,xeps,keyerr,anyerr);
  fprintf(stderr,"xent: %15.6E xcs2: %15.6E xdedt: %15.6E\n",xent,xcs2,xdedt);
  fprintf(stderr,"xdpderho: %15.6E xdpdrhoe: %15.6E \n",xdpderho,xdpdrhoe);
  fprintf(stderr,"xxa: %15.6E xxh: %15.6E xxn: %15.6E xxp: %15.6e\n",
	  xxa,xxh,xxn,xxp);
  fprintf(stderr,"xabar: %15.6E xzbar: %15.6E \n",xabar,xzbar);
  fprintf(stderr,"xmue: %15.6E xmun: %15.6E  xmup: %15.6E  xmuhat: %15.6E\n",
	  xmue,xmun,xmup,xmuhat);


  
  // set up the files
  const double k_MeV = 8.6173324e-11;   // boltzmann constant (Mev/K)
  std::ifstream infile;
  infile.open(argv[1]);
  std::ofstream outfile;
  outfile.open("equilibrium_ye.dat");
  int nlines;
  infile >> nlines;
  double x;

  for(int i=0; i<nlines; i++){
    infile >> x;
    infile >> xrho;
    xrho *= RHOGF;
    infile >> xtemp;
    xtemp *= k_MeV;

    outfile << x << "\t" << beta_eq_ye(xrho,xtemp) << std::endl;
  }

  return 0;
}
