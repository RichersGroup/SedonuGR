#include <cstdio>
#include <cstdlib>
#include <sstream>
#include <string>
#include "nuc_eos.hh"
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

using namespace nuc_eos;
using namespace std;
const double k_MeV = 8.6173324e-11;   // boltzmann constant (Mev/K)
const double mn = 939.56563; // MeV
const double mp = 938.27231; // MeV

int main(int argc, char* argv[]) {

	if(argc!=3){
		std::cout << "Usage: many_Ye EOSfile neutron_star.mod" << std::endl;
		exit(1);
	}

	//read in parameters
	nuc_eos_C_ReadTable(argv[1]);
	ifstream infile(argv[2]);

	// get file properties from first line
	string trash;
	infile >> trash; // grid type
	if(trash != "1D_sphere"){
	  cout << "Error: need 1D_sphere grid." << endl;
	}
	int n_zones;
	infile >> n_zones;
	double rmin;
	infile >> rmin; // rmin
	
	// read the file
	vector<double> rtop(n_zones);
	vector<double> rmid(n_zones);
	vector<double> rho(n_zones);
	vector<double> Ye(n_zones);
	vector<double> T(n_zones);
	for(int z_ind=0; z_ind<n_zones; z_ind++){
	  infile >> rtop[z_ind];
	  infile >> rho[z_ind];
	  infile >> T[z_ind];
	  infile >> Ye[z_ind];
	  infile >> trash; //tmp_vr[z_ind];
	  infile >> trash; //tmp_alpha[z_ind];
	  infile >> trash; //tmp_X[z_ind];
	  rho[z_ind] *= RHOGF; // convert to code units
	  T[z_ind] *= k_MeV; // convert to code units
	  
	  double last = z_ind==0 ? rmin : rtop[z_ind-1];
	  rmid[z_ind] = 0.5 * (rtop[z_ind] + last);

	  
	  // set up the solver
	  double solve_for_munu = 0;//(mn-mp); // eos does not include baryon masses
	  double ye_min = 0.05;
	  double ye_max = 0.55;
	  double xye = 0.0;
	  double xeps, xprs, xent, xcs2, xdedt, xdpderho, xdpdrhoe;
	  double xa, xh, xn, xp, abar, zbar;
	  double mue, mun, mup, muhat, munu=0;
	  int keyerr, anyerr;
	  int n=1;
	  int iter=0;
	  int maxiter = 1000;
	  do{
	    xye = 0.5*(ye_min + ye_max);
	    nuc_eos_m_kt1_full(&n,&rho[z_ind],&T[z_ind],&xye,
			       &xeps,&xprs,&xent,&xcs2,&xdedt,
			       &xdpderho,&xdpdrhoe,&xa,&xh,&xn,&xp,&abar,&zbar,
			       &mue,&mun,&mup,&muhat,
			       &keyerr,&anyerr);
	    munu = mue-muhat;
	    if(munu-solve_for_munu<0) ye_min = xye;
	    else ye_max = xye;
	    iter++;
	  } while(fabs(munu-solve_for_munu) > 1e-10 && iter<maxiter);
	  std::cout << rmid[z_ind] << "\t" << rho[z_ind]/RHOGF << "\t" << T[z_ind]/k_MeV << "\t" << xye << std::endl;
	}

	return 0;
}
