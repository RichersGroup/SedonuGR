#ifndef _ZONE_H
#define _ZONE_H

// define real to choose either double or float precision
typedef double real;


//-------------------------------------------------
// Class to store properties of one zone
//-------------------------------------------------
class zone
{

public:

  // fluid properties
  real v[3];            // velocity vector (cm/s)
  real rho;             // density (g/cm^3)
  real cs;              // sound speed (cm/s)
  real p_gas;           // gas pressure
  real e_gas;           // gas energy density per gram
  real E_gas;           // gas total energy
  real T_gas;           // gas temperature
  real Ye;              // electron fraction

  // store other parameters
  real ni56;            // nickel fraction
  real H;               // specific heating rate (erg/s/g)

  // radiation quantities
  real e_rad;      // radiation energy density  (ergs/cm^3) in lab frame
  real e_abs;      // radiation energy deposition density rate (ergs/cm^3/s)
  real fx_rad;     // radiation x-force in lab frame
  real fy_rad;     // radiation y-force in lab frame
  real fz_rad;     // radiation z-force in lab frame
  real eps_imc;    // fleck factor effective absorption

  real l_abs;      // lepton number deposition density rate (cm^-3 s^-1)

  // four force vector in lab frame
  real G1, G2, G3;

  // radiation pessure tensor components (symmetric)
  real P11, P12, P13, P22, P23, P33;
};

#endif
