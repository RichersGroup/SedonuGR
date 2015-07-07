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

#include "zone.h"
#include <iomanip>
#include <cmath>
#include <fstream>
#include "global_options.h"

zone::zone(const int dimensionality){
	v.resize(dimensionality);
	rho = NaN;
	T = NaN;
	Ye = NaN;
	H_vis = NaN;
	e_abs = NaN;
	nue_abs = NaN;
	anue_abs = NaN;
	e_emit = NaN;
	l_emit = NaN;
	t_eabs = NaN;
	t_eemit = NaN;
	t_labs = NaN;
	t_lemit = NaN;
	Q_annihil = NaN;
}
