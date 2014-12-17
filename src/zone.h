#ifndef _ZONE_H
#define _ZONE_H
#include "global_options.h"
#include <vector>
#include <fstream>
#include "spectrum_array.h"

//-------------------------------------------------
// Class to store properties of one zone
//-------------------------------------------------
class zone
{

public:

	// constructors
	zone(const int dimensionality=0);

	// fluid properties (rho,T are in comoving frame. Ye is invariant.)
	vector<double> v;       // velocity vector (cm/s)
	double rho;             // density (g/cm^3)
	double T;               // gas temperature (K)
	double Ye;              // electron fraction

	// store other parameters
	double H_vis;               // specific heating rate (erg/s/g)

	// radiation quantities (all in comoving frame) (dVdt and lepton number are relativistic invariants)
	double e_rad;                         // radiation energy density  (ergs/cm^3) in lab frame
	double e_abs;                         // radiation energy deposition density rate (ergs/cm^3/s)
	double nue_abs;                       // electron neutrino number deposition density rate (cm^-3 s^-1)
	double anue_abs;                      // electron anti-neutrino number deposition density rate (cm^-3 s^-1)
	double e_emit;                        // radiation energy emission rate (erg/ccm/s)
	double l_emit;						  // lepton number emission rate (cm^-3 s^-1)
	vector<spectrum_array> distribution;  // radiation energy density for each species (erg/ccm. Integrated over bin frequency and direction)
	double Q_annihil;                     // annihilation energy deposition rate (erg/ccm/s)
	double nu_avg;

	// timescales (all in lab frame)
	double t_eabs;    // heating timescale
	double t_eemit;    // cooling timescale
	double t_labs;    // leptonization timescale from absorption
	double t_lemit;   // leptonization timescale from emission
};

#endif

// NOTES
// - tried implementing as struct of vector rather than vector of structs (which would eliminate the need for the mpi datatype definitions). Significantly slower, presumably b/c access patterns are not sequential. Also, this probably increased false sharing, making OpenMP less efficient (guess)
