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

#ifndef _GRID_1D_SPHERE_H
#define _GRID_1D_SPHERE_H 1

#include "Grid.h"

using std::vector;

//*******************************************
// 1-Dimensional Spherical geometry
//*******************************************
class Grid1DSphere: public Grid
{

protected:
	// store location of the outer edge of the zone.
	LocateArray r_out;
	int reflect_outer;

	// ds^2 = -alpha dt^2 + X^2 dr^2 + dOmega^2
	class Metric{
	private:
		int array_len = 0;
		vector<double> alpha;
		vector<double> X;
		vector<double> dadr;
		vector<double> dXdr;

	public:

		int size() const {
			PRINT_ASSERT(alpha.size(),==,array_len);
			PRINT_ASSERT(X.size()    ,==,array_len);
			PRINT_ASSERT(dadr.size() ,==,array_len);
			PRINT_ASSERT(dXdr.size() ,==,array_len);
			return array_len;
		}
		void resize(const int len){
			alpha.resize(len);
			X.resize(len);
			dadr.resize(len);
			dXdr.resize(len);
			array_len = len;
		}
		void set_metric(const vector<double>& tmp_alpha, const vector<double>& tmp_X,
				const vector<double>& tmp_dadr, const vector<double>& tmp_dXdr){
			PRINT_ASSERT(tmp_alpha.size(),==,array_len);
			PRINT_ASSERT(tmp_X.size()    ,==,array_len);
			PRINT_ASSERT(tmp_dadr.size() ,==,array_len);
			PRINT_ASSERT(tmp_dXdr.size() ,==,array_len);

			this->alpha = tmp_alpha;
			this->X     = tmp_X;
			this->dadr  = tmp_dadr;
			this->dXdr  = tmp_dXdr;
		}
		double get_alpha(const int z_ind, const double r, const LocateArray& r_out)const {
			return alpha[z_ind] + dadr[z_ind] * (r-r_out.center(z_ind));}
		double get_X(    const int z_ind, const double r, const LocateArray& r_out) const {
			return     X[z_ind] + dXdr[z_ind] * (r-r_out.center(z_ind));}
		double get_dadr( const int z_ind) const {
			return  dadr[z_ind];}
		double get_dXdr( const int z_ind) const {
			return  dXdr[z_ind];}
	};
	Metric metric;

public:

	Grid1DSphere();
	virtual ~Grid1DSphere() {}

	void read_model_file(Lua* lua);
	void read_custom_model(Lua* lua);
	void read_nagakura_model(Lua* lua);

	// required functions
	int  zone_index               (const double x[3], const int xsize                            ) const;
	double zone_lab_volume        (const int z_ind                                               ) const;
	double zone_min_length        (const int z_ind                                               ) const;
	void zone_coordinates         (const int z_ind, double r[1], const int rsize                 ) const;
	void zone_directional_indices (const int z_ind, int dir_ind[1], const int size               ) const;
	void sample_in_zone (const int z_ind, const double rand[3], const int randsize, double x[3], const int xsize) const;
	void interpolate_fluid_velocity(const double x[3], const int xsize, double v[3], const int vsize, int z_ind) const;
	void write_rays               (const int iw                                                  ) const;
	void symmetry_boundaries      (LorentzHelper *lh                                             ) const;
	double lab_dist_to_boundary   (const LorentzHelper *lh                                       ) const;
	double zone_radius            (const int z_ind) const;
	void dims                     (hsize_t dims[1], const int size) const;
	hsize_t dimensionality() const {return 1;};
	void write_hdf5_coordinates(H5::H5File file) const;
	double zone_cell_dist(const double x_up[3], const int z_ind) const;

	// GR functions
	void g_down(const double xup[4], double g[4][4]) const;
	void connection_coefficients(const double xup[4], double gamma[4][4][4]) const; // Gamma^alhpa_mu_nu
};


#endif
