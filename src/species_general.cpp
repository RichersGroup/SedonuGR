#include "species_general.h"
#include <gsl/gsl_rng.h>
#include <math.h>
#include <cassert>
#include <limits>
#include "transport.h"
#include "physical_constants.h"
#include "Lua.h"
#include <list>

#define NaN std::numeric_limits<double>::quiet_NaN()
#define MAX std::numeric_limits<int>::max()
namespace pc = physical_constants;

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
}

void species_general::init(Lua* lua, transport* simulation)
{
	// set the pointer to see the simulation info
	sim = simulation;

	// call child's init function
	myInit(lua);

	// allocate space for the grid eas spectrum containers
	abs_opac.resize(sim->grid->z.size());
	scat_opac.resize(sim->grid->z.size());
	emis.resize(sim->grid->z.size());

	// allocate space for each eas spectrum
	if(sim->n_emit_core > 0) core_emis.resize(nu_grid.size());
    #pragma omp parallel for
	for (int i=0; i<abs_opac.size();  i++){
		abs_opac[i].resize(nu_grid.size());
		scat_opac[i].resize(nu_grid.size());
		emis[i].resize(nu_grid.size());
	}

	// set up core emission spectrum function (erg/s)
	if(sim->n_emit_core > 0){
		//double L_core = lua->scalar<double>("L_core");
		double T_core = lua->scalar<double>("T_core") / pc::k_MeV;
		double r_core = lua->scalar<double>("r_core");
		double chempot = lua->scalar<double>("core_nue_chem_pot") * (double)lepton_number * pc::MeV_to_ergs;
        #pragma omp parallel for ordered
		for (int j=0;j<nu_grid.size();j++)
		{
			double nu  = nu_grid.center(j);
			double dnu = nu_grid.delta(j);
            #pragma omp ordered
			core_emis.set_value(j, blackbody(T_core,chempot,nu)*dnu);
		}
		core_emis.normalize();
		core_emis.N *= pc::pi * (4.0*pc::pi*r_core*r_core) * weight;
	}
}




//----------------------------------------------------------------
// return a randomly sampled frequency
// for a particle emitted from the core
//----------------------------------------------------------------
double species_general::sample_core_nu() const
{
	// randomly pick a frequency
	double rand = sim->rangen.uniform();
	double result = core_emis.invert_linear(rand,&nu_grid);
	assert(result>0);
	return result;
}

//----------------------------------------------------------------
// return a randomly sampled frequency
// for a particle emitted from a zone
//----------------------------------------------------------------
double species_general::sample_zone_nu(const int zone_index) const
{
	// randomly pick a frequency
	double rand = sim->rangen.uniform();
	double result = emis[zone_index].invert_linear(rand,&nu_grid);
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
	for(int i=0; i<emis[zone_index].size(); i++)
	{
		l_emis += lepton_number * emis[zone_index].get_value(i) / (pc::h*nu_grid.x[i]);
	}
	return l_emis * emis[zone_index].N;
}
