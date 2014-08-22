#include <limits>
#include <vector>
#include <cassert>
#include "photons.h"
#include "transport.h"
#include "Lua.h"
#include "global_options.h"

//----------------------------------------------------------------
// called from species_general::init (photon-specific stuff)
//----------------------------------------------------------------
void photons::myInit(Lua* lua)
{
	// set name
	name = "Photons";
	weight = 1.0;

	// intialize output spectrum
	vector<double>stg = lua->vector<double>("phot_spec_time_grid");
	vector<double>sng = lua->vector<double>("phot_spec_nu_grid");
	int nmu  = lua->scalar<int>("phot_spec_n_mu");
	int nphi = lua->scalar<int>("phot_spec_n_phi");
	spectrum.init(stg,sng,nmu,nphi);

	// read opacity parameters
	grey_opac       = lua->scalar<double>("phot_grey_opacity");
	grey_abs_frac   = lua->scalar<double>("phot_grey_abs_frac");
	double nu_start = lua->scalar<double>("phot_nugrid_start");
	double nu_stop  = lua->scalar<double>("phot_nugrid_stop");
	int      n_nu   = lua->scalar<int>("phot_nugrid_n");
	lepton_number   = 0;
	assert(nu_stop >= nu_start);

	// initialize the  frequency grid
	nu_grid.init(nu_start,nu_stop,n_nu);

	// set photon's min and max values
	T_min   =  1.0;
	T_max   =  1e12;
	Ye_min  = -numeric_limits<double>::infinity();
	Ye_max  =  numeric_limits<double>::infinity();
	rho_min = -numeric_limits<double>::infinity();
	rho_max =  numeric_limits<double>::infinity();
}


//-----------------------------------------------------------------
// set emissivity, abs. opacity, and scat. opacity in zones
//-----------------------------------------------------------------
void photons::set_eas(int zone_index)
{
	zone* z = &(sim->grid->z[zone_index]);
	assert(grey_opac >= 0);

	// leave serial. Parrallelized threads call this function.
	for(unsigned j=0;j<nu_grid.size();j++)
	{
		double nu  = nu_grid.center(j);        // (Hz)
		double dnu = nu_grid.delta(j);         // (Hz)
		double bb  = blackbody(z->T_gas,0,nu)*dnu;  // (erg/s/cm^2/ster)
		emis[zone_index].set_value(j,grey_opac*grey_abs_frac*bb*z->rho); // (erg/s/cm^3/Hz/ster)
	}
	emis[zone_index].normalize();
}



//-----------------------------------------------------------------
// Klein_Nishina correction to the Compton cross-section
// assumes energy x is in MeV
//-----------------------------------------------------------------
double photons::klein_nishina(const double x_input) const
{
	// divide by m_e c^2 = 0.511 MeV
	double x = x_input/pc::m_e_MeV;
	double logfac = log(1 + 2*x);
	double term1 = (1+x)/x/x/x*(2*x*(1+x)/(1+2*x) - logfac);
	double term2 = 1.0/2.0/x*logfac;
	double term3 = -1.0*(1 + 3*x)/(1+2*x)/(1+2*x);
	double KN    = .75*(term1 + term2 + term3);
	return KN;
}

//================================================//
// calculate planck function (erg/s/cm^2/Hz/ster) //
//================================================//
double photons::blackbody(const double T, const double chempot, const double nu) const
{
	assert(chempot==0);
	double zeta = (pc::h*nu - chempot) / (pc::k*T);
	return 2.0*nu*nu*nu*pc::h/pc::c/pc::c/(exp(zeta)-1);
}
