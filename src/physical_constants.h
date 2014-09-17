#ifndef PHYS_CONST
#define PHYS_CONST 1

namespace physical_constants
{
const double pi    = 3.14159;         // just pi
const double c     = 2.99792458e10;   // speed of light (cm/s)
const double h     = 6.6260755e-27;   // plancks constant (ergs s)
const double k     = 1.380658e-16;    // boltz constatn (ergs/K)
const double k_ev  = 8.6173324e-5;    // boltzmann constant (ev/K)
const double m_p   = 1.67262158e-24;  // proton mass (g)
const double m_n   = 1.67492861e-24;  // neutron mass (g)
const double m_e   = 9.10938188e-28;  // mass of electron (g)
const double m_alpha = 6.64465675e-24; // mass of alpha particle (g)
const double sb    = 5.6704e-5;       // stefan boltzman constant (ergs cm^-2 s^-1 K^-4)
const double a     = 7.5657e-15;      // radiation constant (ergs cm-3 K-4)
const double h_MeV = 4.135667516e-21; // planck constant (MeV s)
const double k_MeV = 8.6173324e-11;   // boltzmann constant (Mev/K)
const double G     = 6.67259e-8;      // gravitational constant (cm^3/g/s^2)
const double M_sun = 1.989e33;        // solar mass (g)

const double m_e_MeV =  0.510998910511;    // rest energy of electron in Mev
const double ergs_to_MeV = 1.0/(1.60217646e-6);
const double MeV_to_ergs = 1.60217646e-6;
const double ev_to_ergs =  1.60217646e-12;
const double cm_to_angs =  1.0e8;
const double angs_to_cm =  1.0e-8;

const double thomson_cs = 0.66523e-24;       // Thomson cross-section cm-2
const double sigma_tot = 0.0265400193567;    // integrated line coefficent (cm^2 Hz)
}


#endif
