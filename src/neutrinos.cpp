#include "global_options.h"
#include <vector>
#include "neutrinos.h"
#include "transport.h"
#include "Lua.h"
#include "nulib_interface.h"
#include <cmath>

// constructor
neutrinos::neutrinos(){
	num_nut_species = MAX;
	nulibID = MAX;
	cutoff=0;
}


//----------------------------------------------------------------
// called from species_general::init (neutrino-specific stuff)
//----------------------------------------------------------------
void neutrinos::myInit(Lua* lua)
{
	// intialize output spectrum
	std::vector<double>sng = lua->vector<double>("nut_spec_nu_grid");
	int nmu  = lua->scalar<int>("nut_spec_n_mu");
	int nphi = lua->scalar<int>("nut_spec_n_phi");
	spectrum.init(sng,nmu,nphi);

	// set lepton number
	if(nulibID == 0)   lepton_number =  1;
	if(nulibID == 1)   lepton_number = -1;
	if(nulibID >= 2)   lepton_number =  0;

	// set names and spectrum weight
	if     (nulibID == 0) {name = "Electron Neutrinos";              weight = 1.0;}
	else if(nulibID == 1) {name = "Electron Anti-Neutrinos";         weight = 1.0;}
	else if(nulibID == 2){
		if     (num_nut_species == 3) {name = "Mu/Tau Anti/Neutrinos"; weight = 4.0;}
		else if(num_nut_species == 4) {name = "Mu/Tau Neutrinos";      weight = 2.0;}
		else if(num_nut_species == 6) {name = "Mu Neutrinos";          weight = 1.0;}
		else                          {name = "ERROR";                 weight = -1;}}
	else if(nulibID == 3){
		if     (num_nut_species == 4) {name = "Mu/Tau Anti-Neutrinos"; weight = 2.0;}
		else if(num_nut_species == 6) {name = "Mu Antineutrino";       weight = 1.0;}
		else                          {name = "ERROR";                 weight = -1;}}
	else if(nulibID == 4){
		if(num_nut_species == 6)      {name = "Tau Neutrinos";         weight = 1.0;}
		else                          {name = "ERROR";                 weight = -1;}}
	else if(nulibID == 5){
		if(num_nut_species == 6)      {name = "Tau Anti-Neutrinos";    weight = 1.0;}
		else                          {name = "ERROR";                 weight = -1;}}
	else{
		cout << "ERROR: Sedona does not know how to deal with a neutrino ID of " << nulibID << "." << endl;
		exit(16);}

	// set up the frequency table
	cutoff        = lua->scalar<double>("nut_cdf_cutoff");
	grey_opac     = lua->scalar<double>("nut_grey_opacity");
	grey_abs_frac = lua->scalar<double>("nut_grey_abs_frac");
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
		double nu_start = lua->scalar<double>("nut_nugrid_start");
		double nu_stop  = lua->scalar<double>("nut_nugrid_stop");
		assert(nu_stop > nu_start);
		int      n_nu   = lua->scalar<int>("nut_nugrid_n");
		nu_grid.init(nu_start,nu_stop,n_nu);

		// set neutrino's min and max values
		T_min   =  1.0;
		T_max   =  1e12;
		Ye_min  = 0;
		Ye_max  = 1;
		rho_min = -numeric_limits<double>::infinity();
		rho_max =  numeric_limits<double>::infinity();
	}

}


//-----------------------------------------------------------------
// set emissivity, abs. opacity, and scat. opacity in zones
//-----------------------------------------------------------------
void neutrinos::set_eas(int zone_index)
{
	zone* z = &(sim->grid->z[zone_index]);
	double ngroups = (double)emis[zone_index].size();

	if(grey_opac < 0){ // get opacities and emissivity from NuLib
		nulib_get_eas_arrays(z->rho, z->T, z->Ye, nulibID,
				emis[zone_index], abs_opac[zone_index], scat_opac[zone_index]);
		emis[zone_index].normalize(cutoff/ngroups);
	}

	else{ // get emissivity from blackbody and the grey opacity
		assert(grey_abs_frac>=0 && grey_abs_frac<=1.0);
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
		}
		emis[zone_index].normalize(cutoff/ngroups);
	}
}


//-----------------------------------------------------------------
// Calculate the fermi-dirac blackbody function (erg/s/cm^2/Hz/ster)
//-----------------------------------------------------------------
double neutrinos::blackbody(const double T, const double chem_pot, const double nu) const
{
	double zeta = (pc::h*nu - chem_pot)/pc::k/T;
	double bb = (pc::h*nu) * (nu/pc::c) * (nu/pc::c) / (exp(zeta) + 1.0);
	assert(bb >= 0);
	return bb;
}

//-------------------------------------------------------------------------------------------------------------
// Calculate the neutrino anti-neutrino annihilation rate (erg/ccm/s) from the distribution functions (erg/ccm)
//-------------------------------------------------------------------------------------------------------------
double cos_angle_between(const double mu1, const double mu2, const double phi1, const double phi2){
	assert(mu1<=1.0 && mu1>=-1.0);
	assert(mu2<=1.0 && mu2>=-1.0);
	double result = mu1*mu2 + sqrt((1.0-mu1*mu1)*(1.0-mu2*mu2))*cos(phi2-phi1); // Bruenn 1985 eq. A7

	// make sure numerical error doesn't push it beyond {-1,1}
	result = min(result, 1.0);
	result = max(result,-1.0);
	return result;
}
double neutrinos::annihilation_rate(
		const spectrum_array& nu_dist,     // erg/ccm (integrated over angular bin and energy bin)
		const spectrum_array& nubar_dist,  // erg/ccm (integrated over angular bin and energy bin)
		const bool electron_type){         // is this an electron-type interaction?

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

	// integrate over all bins
	double Q = 0;
	for(unsigned i=0; i<nu_dist.size(); i++){
		double avg_e   = nu_dist.nu_center(i)*pc::h; // erg
		double avg_mu  = nu_dist.mu_center(i);
		double avg_phi = nu_dist.phi_center(i);
		double nudist_edens = nu_dist.get(i); // erg/ccm

		for(unsigned j=0; j<nubar_dist.size(); j++){
			double avg_ebar   = nubar_dist.nu_center(j)*pc::h; // erg
			double avg_mubar  = nubar_dist.mu_center(j);
			double avg_phibar = nubar_dist.phi_center(j);
			double nubardist_edens = nubar_dist.get(j); // erg/ccm

			double costheta = cos_angle_between(avg_mu,avg_mubar,avg_phi,avg_phibar);
			Q += nudist_edens * nubardist_edens * (avg_e + avg_ebar) * (1.0-costheta) * (
					C1pC2/3.0 * (1-costheta) +
					C3 * mec2*mec2 / (avg_e*avg_ebar)
					);
		}
	}

	Q *= pc::sigma0*pc::c / (4.*mec2*mec2); // erg/ccm/s
	assert(Q >= 0);

	return Q;
}
