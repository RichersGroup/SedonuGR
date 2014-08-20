#ifndef _ZONE_H
#define _ZONE_H
#include <vector>
#include <fstream>
using namespace std;

// define real to choose either double or float precision
typedef double real;

//-------------------------------------------------
// Class to store properties of one zone
//-------------------------------------------------
class zone
{

public:

  // constructors
  zone(const int dimensionality=0);

  // writing utilities
  static void write_header(const int dimensionality, ofstream& outf);
  void write_line(const vector<double>& r, ofstream& outf) const;

  // fluid properties
  vector<real> v;            // velocity vector (cm/s)
  real rho;             // density (g/cm^3)
  real T_gas;           // gas temperature (K)
  real Ye;              // electron fraction

  // store other parameters
  real H;               // specific heating rate (erg/s/g)

  // radiation quantities
  real e_rad;                         // radiation energy density  (ergs/cm^3) in lab frame
  real e_abs;                         // radiation energy deposition density rate (ergs/cm^3/s)
  real l_abs;                         // lepton number deposition density rate (cm^-3 s^-1)
  real e_emit;                        // radiation energy emission rate (erg/ccm/s)
  real l_emit;						  // lepton number emission rate (cm^-3 s^-1)

  // timescales
  real t_eabs;    // heating timescale
  real t_eemit;    // cooling timescale
  real t_labs;    // leptonization timescale from absorption
  real t_lemit;   // leptonization timescale from emission
};

#endif

// NOTES
// - tried implementing as struct of vector rather than vector of structs (which would eliminate the need for the mpi datatype definitions). Significantly slower, presumably b/c access patterns are not sequential. Also, this probably increased false sharing, making OpenMP less efficient (guess)
