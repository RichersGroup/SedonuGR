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

#include "SpectrumArray.h"
#include "global_options.h"

SpectrumArray::SpectrumArray(){
	rotated_basis = -MAXLIM;
}

void SpectrumArray::rotate_basis(double D[3], const double xyz[3]){
	double x=xyz[0], y=xyz[1], z=xyz[2];
	double inv_r = 1.0 / sqrt(Grid::dot_Minkowski<3>(xyz,xyz));
	double rp = sqrt(x*x + y*y);
	double rhat[3]     = {0,0,0};
	double thetahat[3] = {0,0,0};
	double phihat[3]   = {0,0,0};
	if(rp==0){
		rhat[2] = z>0 ? -1.0 : 1.0;
	    thetahat[1] = 1;
	    phihat[0] = 1;
	}
	else{
		double inv_rp = 1.0/rp;
		rhat[0] = x*inv_r;
		rhat[1] = y*inv_r;
		rhat[2] = z*inv_r;
		thetahat[0] = z*inv_r * x*inv_rp;
		thetahat[1] = z*inv_r * y*inv_rp;
		thetahat[2] = -rp * inv_r;
		phihat[0] = -y*inv_rp;
		phihat[1] =  x*inv_rp;
		phihat[2] = 0;
	}

	double Dout[3];
	Dout[0] = Grid::dot_Minkowski<3>(D,thetahat);
	Dout[1] = Grid::dot_Minkowski<3>(D,phihat);
	Dout[2] = Grid::dot_Minkowski<3>(D,rhat);
	for(int i=0; i<3; i++) D[i] = Dout[i];
}

void SpectrumArray::rotate_and_count(const double kup[4], const double xup[4], const double nu, const double weight){
	double D[3] = {kup[0], kup[1], kup[2]};
	if(rotated_basis) rotate_basis(D,xup);

	Grid::normalize_Minkowski<3>(D);
	count(D, nu, weight);
}
