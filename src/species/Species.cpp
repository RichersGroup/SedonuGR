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

#include <mpi.h>
#include <fstream>
#include "global_options.h"
#include "physical_constants.h"
#include "Species.h"
#include "PolarSpectrumArray.h"
#include "MomentSpectrumArray.h"
#include "RadialMomentSpectrumArray.h"
#include "GR1DSpectrumArray.h"
#include "Transport.h"
#include "Grid.h"

using namespace std;
namespace pc = physical_constants;

Species::Species(){
	weight = NaN;
	lepton_number = MAXLIM;
	num_species = MAXLIM;
	ID = MAXLIM;
	T_core = NaN;
	mu_core = NaN;
	core_lum_multiplier = NaN;
}

void Species::init(Lua* lua)
{
	// set lepton number
	if(ID == 0)   lepton_number =  1;
	if(ID == 1)   lepton_number = -1;
	if(ID >= 2)   lepton_number =  0;

	// set names and spectrum weight
	if     (ID == 0) {name = "Electron Neutrinos";              weight = 1.0;}
	else if(ID == 1) {name = "Electron Anti-Neutrinos";         weight = 1.0;}
	else if(ID == 2){
		if     (num_species == 3) {name = "Mu/Tau Anti/Neutrinos"; weight = 4.0;}
		else if(num_species == 4) {name = "Mu/Tau Neutrinos";      weight = 2.0;}
		else if(num_species == 6) {name = "Mu Neutrinos";          weight = 1.0;}
		else                          {name = "ERROR";                 weight = -1;}}
	else if(ID == 3){
		if     (num_species == 4) {name = "Mu/Tau Anti-Neutrinos"; weight = 2.0;}
		else if(num_species == 6) {name = "Mu Antineutrino";       weight = 1.0;}
		else                          {name = "ERROR";                 weight = -1;}}
	else if(ID == 4){
		if(num_species == 6)      {name = "Tau Neutrinos";         weight = 1.0;}
		else                          {name = "ERROR";                 weight = -1;}}
	else if(ID == 5){
		if(num_species == 6)      {name = "Tau Anti-Neutrinos";    weight = 1.0;}
		else                          {name = "ERROR";                 weight = -1;}}
	else{
		cout << "ERROR: Sedona does not know how to deal with a neutrino ID of " << ID << "." << endl;
		exit(16);}

	//============================//
	// CALL CHILD'S INIT FUNCTION //
	//============================//
	myInit(lua);
}

// cm^3/s
void Species::get_annihil_kernels(const double /*rho*/, const double /*T*/, const double /*Ye*/, const Axis& nuAxis, vector< vector< vector<double> > >& phi) const{
	// constants
	using namespace pc;
	double mec2 = m_e*c*c; // erg
	double C1pC2,C3; // from Ruffert+97 (paper II) above equation 4
	if(ID==0 || ID==1){ // electron-type
		C1pC2 = (CV-CA)*(CV-CA) + (CV+CA)*(CV+CA); // ~2.34
		C3 = 2./3. * (2.*CV*CV - CA*CA);           // ~1.06
	}
	else{
		C1pC2 = (CV-CA)*(CV-CA) + (CV+CA-2.)*(CV+CA-2.);     // ~0.50
		C3 = 2./3. * (2.*(CV-1.)*(CV-1.) - (CA-1.)*(CA-1.)); // ~-0.16
	}
	double twomec2 = 2.0*pc::m_e*pc::c*pc::c;
	double C3mec4 = C3*mec2*mec2;
	double C1pC2_3 = C1pC2/3.0;
	double mec22 = mec2*mec2;
	size_t nnu = nuAxis.size();

	phi.resize(3);
	for(size_t k=0; k<3; k++) phi[k].resize(nnu);

	for(size_t inu=0; inu<nnu; inu++){
		double avg_e = nuAxis.mid[inu]*pc::h; // erg
		for(size_t k=0; k<3; k++) phi[k][inu].resize(nnu);

		for(size_t inubar=0; inubar<nnu; inubar++){
			double avg_ebar = nuAxis.mid[inubar]*pc::h; // erg

			double eebar = avg_e * avg_ebar;
			double C3mec4_eebar = C3mec4/eebar;
			// phi(mu) = A( B(1-mu)^2 + C(1-mu))
			double A = pc::sigma0*pc::c * eebar/(twomec2*twomec2); // cm^3/s
			double B = C1pC2_3;
			double C = C3mec4_eebar;
			if(eebar > mec22){
				phi[0][inu][inubar] =  2.*A*(4.*B/3. + C);
				phi[1][inu][inubar] = -2./3.*A*(2.*B + C);
				phi[2][inu][inubar] =  4./15.*A*B;
				PRINT_ASSERT(phi[0][inu][inubar],>=,0);
				// sanity checks. format: coeff*phi*P(n,x)
				PRINT_ASSERT((1./2.*phi[0][inu][inubar] + 3./2.*phi[1][inu][inubar]*( 1  ) + 5./2.*phi[2][inu][inubar]*(1   ))/phi[0][inu][inubar],>=,-TINY); // x= 1
				PRINT_ASSERT((1./2.*phi[0][inu][inubar] + 3./2.*phi[1][inu][inubar]*(-1  ) + 5./2.*phi[2][inu][inubar]*(1   ))/phi[0][inu][inubar],>=,-TINY); // x=-1
				PRINT_ASSERT((1./2.*phi[0][inu][inubar] + 3./2.*phi[1][inu][inubar]*( 0  ) + 5./2.*phi[2][inu][inubar]*(-0.5))/phi[0][inu][inubar],>=,-TINY); // x= 0

			}
		}
	}
}
