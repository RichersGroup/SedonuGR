/*
//  Copyright (c) 2015, California Institute of Technology and the Regents
//  of the University of California, based on research sponsored by the
//  United States Department of Energy. All rights reserved.
//
//  This file is part of Sedonu.
//
//  Sedonu is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  Neither the name of the California Institute of Technology (Caltech)
//  nor the University of California nor the names of its contributors 
//  may be used to endorse or promote products derived from this software
//  without specific prior written permission.
//
//  Sedonu is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with Sedonu.  If not, see <http://www.gnu.org/licenses/>.
//
*/

#include "global_options.h"
#include "Neutrino.h"
#include "Transport.h"
#include "Grid.h"
#include "nulib_interface.h"

using namespace std;
namespace pc = physical_constants;

// constructor
Neutrino::Neutrino(){
	num_species = MAXLIM;
	nulibID = MAXLIM;
	cutoff=0;
}


//----------------------------------------------------------------
// called from species_general::init (neutrino-specific stuff)
//----------------------------------------------------------------
void Neutrino::myInit(Lua* lua)
{
	// set up the frequency table
	cutoff        = lua->scalar<double>("cdf_cutoff");
	grey_opac     = lua->scalar<double>("grey_opacity");
	if(grey_opac < 0){
		nulib_get_nu_grid(nu_grid);

		// set neutrino's min and max values
		T_min  =  nulib_get_Tmin();
		T_max  =  nulib_get_Tmax();
		Ye_min =  nulib_get_Yemin();
		Ye_max =  nulib_get_Yemax();
		rho_min = nulib_get_rhomin();
		rho_max = nulib_get_rhomax();

	}
	else{
		grey_abs_frac = lua->scalar<double>("grey_abs_frac");
		double nu_start = lua->scalar<double>("nugrid_start");
		double nu_stop  = lua->scalar<double>("nugrid_stop");
		PRINT_ASSERT(nu_stop,>,nu_start);
		int      n_nu   = lua->scalar<int>("nugrid_n");
		nu_grid.init(nu_start/pc::h_MeV,nu_stop/pc::h_MeV,n_nu);

		// set neutrino's min and max values
		T_min   =  0.0;
		T_max   =  numeric_limits<double>::infinity();
		Ye_min  = 0;
		Ye_max  = 1;
		rho_min = -numeric_limits<double>::infinity();
		rho_max =  numeric_limits<double>::infinity();
	}

	// intialize output spectrum
	PRINT_ASSERT(nu_grid.size(),>,0);
	int nmu  = lua->scalar<int>("spec_n_mu");
	int nphi = lua->scalar<int>("spec_n_phi");
	SpectrumArray tmp_spectrum;
	LocateArray tmp_mugrid, tmp_phigrid;
	tmp_mugrid.init( -1     , 1     , nmu );
	tmp_phigrid.init(-pc::pi, pc::pi, nphi);
	spectrum.init(nu_grid, tmp_mugrid, tmp_phigrid);

	// set lepton number
	if(nulibID == 0)   lepton_number =  1;
	if(nulibID == 1)   lepton_number = -1;
	if(nulibID >= 2)   lepton_number =  0;

	// set names and spectrum weight
	if     (nulibID == 0) {name = "Electron Neutrinos";              weight = 1.0;}
	else if(nulibID == 1) {name = "Electron Anti-Neutrinos";         weight = 1.0;}
	else if(nulibID == 2){
		if     (num_species == 3) {name = "Mu/Tau Anti/Neutrinos"; weight = 4.0;}
		else if(num_species == 4) {name = "Mu/Tau Neutrinos";      weight = 2.0;}
		else if(num_species == 6) {name = "Mu Neutrinos";          weight = 1.0;}
		else                          {name = "ERROR";                 weight = -1;}}
	else if(nulibID == 3){
		if     (num_species == 4) {name = "Mu/Tau Anti-Neutrinos"; weight = 2.0;}
		else if(num_species == 6) {name = "Mu Antineutrino";       weight = 1.0;}
		else                          {name = "ERROR";                 weight = -1;}}
	else if(nulibID == 4){
		if(num_species == 6)      {name = "Tau Neutrinos";         weight = 1.0;}
		else                          {name = "ERROR";                 weight = -1;}}
	else if(nulibID == 5){
		if(num_species == 6)      {name = "Tau Anti-Neutrinos";    weight = 1.0;}
		else                          {name = "ERROR";                 weight = -1;}}
	else{
		cout << "ERROR: Sedona does not know how to deal with a neutrino ID of " << nulibID << "." << endl;
		exit(16);}



}


//-----------------------------------------------------------------
// set emissivity, abs. opacity, and scat. opacity in zones
//-----------------------------------------------------------------
void Neutrino::set_eas(int zone_index)
{
	Zone* z = &(sim->grid->z[zone_index]);
	double ngroups = (double)emis[zone_index].size();

	if(grey_opac < 0){ // get opacities and emissivity from NuLib
		nulib_get_eas_arrays(z->rho, z->T, z->Ye, nulibID,
				emis[zone_index], abs_opac[zone_index], scat_opac[zone_index]);

		// set the biased emissivity
		for(int g=0; g<nu_grid.size(); g++)
			biased_emis[zone_index].set_value(g, emis[zone_index].get_value(g)
				* sim->importance(abs_opac[zone_index][g], scat_opac[zone_index][g], sim->grid->zone_min_length(zone_index)));
	}

	else{ // get emissivity from blackbody and the grey opacity
		PRINT_ASSERT(grey_abs_frac,>=,0);
		PRINT_ASSERT(grey_abs_frac,<=,1.0);
		for(unsigned j=0;j<nu_grid.size();j++)
		{
			double nu  = nu_grid.center(j);        // (Hz)
			double dnu = nu_grid.delta(j);         // (Hz)
			double bb  = blackbody(z->T,0*pc::MeV_to_ergs,nu)*dnu;  // (erg/s/cm^2/ster)

			double a = grey_opac*z->rho*grey_abs_frac;
			double s = grey_opac*z->rho*(1.0-grey_abs_frac);

			emis[zone_index].set_value(j,a*bb); // (erg/s/cm^3/ster)
			abs_opac[zone_index][j] = a;        // (1/cm)
			scat_opac[zone_index][j] = s;       // (1/cm)

			// set the biased emissivity
			biased_emis[zone_index].set_value(j, emis[zone_index].get_value(j)
					* sim->importance(a, s, sim->grid->zone_min_length(zone_index)));
		}
	}

	emis[zone_index].normalize(cutoff/(double)ngroups);
	biased_emis[zone_index].normalize(cutoff/(double)ngroups);
}


//-----------------------------------------------------------------
// Calculate the fermi-dirac blackbody function (erg/s/cm^2/Hz/ster)
//-----------------------------------------------------------------
double Neutrino::blackbody(const double T /*K*/, const double chem_pot /*erg*/, const double nu /*Hz*/) const
{
	double zeta = (pc::h*nu - chem_pot)/pc::k/T;
	double nu_c = nu*pc::inv_c;
	double bb = (pc::h*nu) * nu_c * nu_c / (exp(zeta) + 1.0);
	PRINT_ASSERT(bb,>=,0);
	return bb;
}

//-------------------------------------------------------------------------------------------------------------
// Calculate the neutrino anti-neutrino annihilation rate (erg/ccm/s) from the distribution functions (erg/ccm)
//-------------------------------------------------------------------------------------------------------------
double cos_angle_between(const double mu1, const double mu2, const double phi1, const double phi2){
	PRINT_ASSERT(mu1,<=,1.0);
	PRINT_ASSERT(mu1,>=,-1.0);
	PRINT_ASSERT(mu2,<=,1.0);
	PRINT_ASSERT(mu2,>=,-1.0);
	double result = mu1*mu2 + sqrt((1.0-mu1*mu1)*(1.0-mu2*mu2))*cos(phi2-phi1); // Bruenn 1985 eq. A7

	// make sure numerical error doesn't push it beyond {-1,1}
	result = min(result, 1.0);
	result = max(result,-1.0);
	return result;
}
double Neutrino::annihilation_rate(
		const SpectrumArray& nu_dist,     // erg/ccm (integrated over angular bin and energy bin)
		const SpectrumArray& nubar_dist,  // erg/ccm (integrated over angular bin and energy bin)
		const bool electron_type,		   // is this an electron-type interaction?
		const int weight){         		   // weight of each species

	PRINT_ASSERT(nu_dist.size(),==,nubar_dist.size());

	// constants
	using namespace pc;
	double mec2 = m_e*c*c; // erg
	double C1pC2,C3; // from Ruffert+97 (paper II) above equation 4
	if(electron_type){
		C1pC2 = (CV-CA)*(CV-CA) + (CV+CA)*(CV+CA); // ~2.34
		C3 = 2./3. * (2.*CV*CV - CA*CA);           // ~1.06
	}
	else{
		C1pC2 = (CV-CA)*(CV-CA) + (CV+CA-2.)*(CV+CA-2.);     // ~0.50
		C3 = 2./3. * (2.*(CV-1.)*(CV-1.) - (CA-1.)*(CA-1.)); // ~-0.16
	}

	// calculate angle between distribution function angles beforehand
	double onemcostheta [nu_dist.mu_dim()] [nu_dist.phi_dim()] [nubar_dist.mu_dim()] [nubar_dist.phi_dim()];
	for(unsigned mu=0; mu<nu_dist.mu_dim(); mu++){
		for(unsigned phi=0; phi<nu_dist.phi_dim(); phi++){

			unsigned index = nu_dist.index(0,mu,phi);
			double avg_mu  = nu_dist.mu_center(index);
			double avg_phi = nu_dist.phi_center(index);

			for(unsigned mubar=0; mubar<nu_dist.mu_dim(); mubar++){
				for(unsigned phibar=0; phibar<nu_dist.phi_dim(); phibar++){

					unsigned indexbar = nubar_dist.index(0,mubar,phibar);
					double avg_mubar  = nubar_dist.mu_center(indexbar);
					double avg_phibar = nubar_dist.phi_center(indexbar);

					onemcostheta[mu][phi][mubar][phibar] = 1.0-cos_angle_between(avg_mu,avg_mubar,avg_phi,avg_phibar);
				}
			}
		}
	}

	// some useful constants
	double twomec2 = 2.0*pc::m_e*pc::c*pc::c;
	double C3mec4 = C3*mec2*mec2;
	double C1pC2_3 = C1pC2/3.0;
	double mec22 = pc::m_e*pc::m_e * pc::c*pc::c*pc::c*pc::c;
	double inv_weight = 1./(double)weight;
	unsigned nnu=nu_dist.nu_dim(), nmu=nu_dist.mu_dim(), nphi=nu_dist.phi_dim();
	unsigned nnubar=nubar_dist.nu_dim(), nmubar=nubar_dist.mu_dim(), nphibar=nubar_dist.phi_dim();

	// integrate over all bins
	double Q = 0;
	// energy loops
	for(unsigned inu=0; inu<nnu; inu++){
		double avg_e = nu_dist.nu_bin_center(inu)*pc::h; // erg
		for(unsigned inubar=0; inubar<nnubar; inubar++){
			double avg_ebar = nubar_dist.nu_bin_center(inubar)*pc::h; // erg

			double C3mec4_eebar = C3mec4/(avg_e*avg_ebar);
			double sume_m_twomec2 = avg_e + avg_ebar - twomec2;
			double eebar = avg_e * avg_ebar;
			if(avg_e*avg_ebar > mec22){

				// neutrino direction loops
				for(unsigned imu=0; imu<nmu; imu++){
					for(unsigned iphi=0; iphi<nphi; iphi++){
						unsigned index = nu_dist.index(inu,imu,iphi);

						// antineutrino direction loops
						for(unsigned imubar=0; imubar<nmubar; imubar++){
							for(unsigned iphibar=0; iphibar<nphibar; iphibar++){
								unsigned indexbar = nubar_dist.index(inubar,imubar,iphibar);

								double onemcost = onemcostheta[imu][iphi][imubar][iphibar];
								if(eebar*onemcost > mec22){
									double nudist_edens    =    nu_dist.get(index)    * inv_weight; // erg/ccm
									double nubardist_edens = nubar_dist.get(indexbar) * inv_weight; // erg/ccm

									Q += nudist_edens * nubardist_edens * sume_m_twomec2 * onemcost *
											(C1pC2_3 * onemcost + C3mec4_eebar);
								} // if
							} // phibar
						} // mubar
					} // phi
				} // mu
			} // if
		} // nubar
	} // nu


	Q *= pc::sigma0*pc::c / (4.*mec2*mec2); // erg/ccm/s
	PRINT_ASSERT(Q,>=,0);

	return Q;
}
