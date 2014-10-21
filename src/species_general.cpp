#include "global_options.h"
#include "species_general.h"
#include <gsl/gsl_rng.h>
#include <math.h>
#include "transport.h"
#include "Lua.h"
#include <list>

species_general::species_general(){
	weight = NaN;
	grey_opac = NaN;
	grey_abs_frac = NaN;
	lepton_number = MAX;
	T_min = NaN;
	T_max = NaN;
	Ye_min = NaN;
	Ye_max = NaN;
	rho_min = NaN;
	rho_max = NaN;
	sim = NULL;
	cdf_interpolation_order = 0;
}

void species_general::init(Lua* lua, transport* simulation)
{
	// set the pointer to see the simulation info
	sim = simulation;
	assert(sim->grid->z.size()>0);

	// call child's init function
	myInit(lua);

	// allocate space for the grid eas spectrum containers
	abs_opac.resize(sim->grid->z.size());
	scat_opac.resize(sim->grid->z.size());
	emis.resize(sim->grid->z.size());

	// allocate space for each eas spectrum
	if(sim->n_emit_core > 0) core_emis.resize(nu_grid.size());
    #pragma omp parallel for
	for(unsigned i=0; i<abs_opac.size();  i++){
		abs_opac[i].resize(nu_grid.size());
		scat_opac[i].resize(nu_grid.size());
		emis[i].resize(nu_grid.size());
	}

	// set up core emission spectrum function (erg/s)
	//double L_core = lua->scalar<double>("L_core");
	double T_core = lua->scalar<double>("T_core") / pc::k_MeV;
	double r_core = lua->scalar<double>("r_core");
	double chempot = lua->scalar<double>("core_nue_chem_pot") * (double)lepton_number * pc::MeV_to_ergs;
	if(sim->n_emit_core > 0) set_cdf_to_BB(T_core, chempot, core_emis);
	core_emis.N *= pc::pi * (4.0*pc::pi*r_core*r_core) * weight;
	cdf_interpolation_order = lua->scalar<int>("cdf_interpolation_order");

    // set up the spectrum in each zone
    bool do_distribution = lua->scalar<int>("do_distribution");
    if(do_distribution){
    	int n_mu = lua->scalar<int>("distribution_nmu");
    	int n_phi = lua->scalar<int>("distribution_nphi");
		#pragma omp parallel for
    	for(unsigned z_ind=0; z_ind<sim->grid->z.size(); z_ind++){
    		spectrum_array tmp_spectrum;
    		locate_array tmp_mugrid, tmp_phigrid;
    		tmp_mugrid.init(-1,1,n_mu);
    		tmp_phigrid.init(0,2.0*pc::pi,n_phi);
    		tmp_spectrum.init(nu_grid, tmp_mugrid, tmp_phigrid);
    		sim->grid->z[z_ind].distribution.push_back(tmp_spectrum);
    	}
    }
}


//-----------------------------------------------------
// set cdf to blackbody distribution
// units of emis.N: erg/s/cm^2/ster
//-----------------------------------------------------
void species_general::set_cdf_to_BB(const double T, const double chempot, cdf_array& emis){
    #pragma omp parallel for ordered
	for(unsigned j=0;j<nu_grid.size();j++)
	{
		double nu  = nu_grid.center(j);
		double dnu = nu_grid.delta(j);
        #pragma omp ordered
		emis.set_value(j, blackbody(T,chempot,nu)*dnu);
	}
	emis.normalize();
}


void species_general::set_emis_to_BB_edens(const double T, const double chempot){
	// This sort of abuses the original meaning of the emissivity arrays for initial particle creation.
	// usage here: units of emis.N are erg/ccm
	// original usage: units of emis.N are erg/s
	for(unsigned z_ind=0; z_ind<emis.size(); z_ind++){
		set_cdf_to_BB(T,chempot,emis[z_ind]);
		emis[z_ind].N *= 4.0*pc::pi/pc::c;
	}
}

//----------------------------------------------------------------
// return a randomly sampled frequency
// for a particle emitted from the core
//----------------------------------------------------------------
double species_general::sample_core_nu() const
{
  assert(nu_grid.min >= 0);

	// randomly pick a frequency
	double rand = sim->rangen.uniform();

	// make sure we don't get a zero frequency
	if(rand==0 && nu_grid.min==0) while(rand==0) rand = sim->rangen.uniform();

	double result = 0;
	assert(cdf_interpolation_order==1 || cdf_interpolation_order==3);
	if     (cdf_interpolation_order==1) result = core_emis.invert_linear(rand,&nu_grid);
	else if(cdf_interpolation_order==3) result = core_emis.invert_cubic(rand,&nu_grid);
	assert(result>0);
	return result;
}

//----------------------------------------------------------------
// return a randomly sampled frequency
// for a particle emitted from a zone
//----------------------------------------------------------------
double species_general::sample_zone_nu(const int zone_index) const
{
  assert(nu_grid.min >= 0);

	// randomly pick a frequency
	double rand = sim->rangen.uniform();

	// make sure we don't get a zero frequency
	if(rand==0 && nu_grid.min==0) while(rand==0) rand = sim->rangen.uniform();

	double result = 0;
	assert(cdf_interpolation_order==1 || cdf_interpolation_order==3);
	if     (cdf_interpolation_order==1) result = emis[zone_index].invert_linear(rand,&nu_grid);
	else if(cdf_interpolation_order==3) result = emis[zone_index].invert_cubic(rand,&nu_grid);
	assert(result>0);
	return result;
}


//----------------------------------------------------------------
// return the emissivity integrated over nu for the core (erg/s)
//----------------------------------------------------------------
double species_general::integrate_core_emis() const
{
	return core_emis.N;
}

//----------------------------------------------------------------
// return the emissivity integrated over nu for a zone (erg/s/ster/cm^3)
//----------------------------------------------------------------
double species_general::integrate_zone_emis(const int zone_index) const
{
	return emis[zone_index].N;
}


//----------------------------------------------------------------
// return the lepton emissivity integrated over nu for a zone (#/s/ster/cm^3)
// ASSUMES linear cdf sampling
//----------------------------------------------------------------
double species_general::integrate_zone_lepton_emis(const int zone_index) const
{
	double l_emis = 0;
	for(unsigned i=0; i<emis[zone_index].size(); i++)
	{
		l_emis += lepton_number * emis[zone_index].get_value(i) / (pc::h*nu_grid.x[i]);
	}
	return l_emis * emis[zone_index].N;
}
