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

//------------------------------------------------------------------
//*****************************************************************
//*************************  GRID ********************************
//*****************************************************************
// The grid class is a construct whose main purpose is to handle 
// the geometry of the model.  It does two main things: (1) Reads
// in the input density,temperature,composition files (that of 
// course must have a specific geometry). (2) Given a set of 
// 3-d coordinates, it will give the corosponding zone index (or
// note that the coords are off the grid).
//
// The grid holds an array of zones, where key fluid data is stored
//
// The grid class is an abstract class that will be used to
// create subclasses (e.g. grid_1D_sphere, grid_3D_cart, etc...
//*****************************************************************


#ifndef _GRID_GENERAL_H
#define _GRID_GENERAL_H 1

#include "Lua.h"
#include "Particle.h"
#include "ThreadRNG.h"
#include "H5Cpp.h"
#include "Axis.h"
#include "MultiDArray.h"
#include "SpectrumArray.h"
#include "Metric.h"
#include "EinsteinHelper.h"

class Transport;
class Zone;
class SpectrumArray;

class Grid
{

protected:

	// fill the grid with data from a model file
	virtual void read_model_file(Lua* lua) = 0;

	// output options
	int do_annihilation;
	int do_GR;

	// get the coordinates at the center of the zone z_ind (GRID COORDINATES)
	virtual void zone_coordinates(const int z_ind, double r[], const int rsize) const = 0;

	// GR functions
	virtual void connection_coefficients(const double xup[4], double gamma[4][4][4], int z_ind=-1) const = 0; // Gamma^alhpa_mu_nu

	double dot (const double aup[4], const double bup[4], const double xup[4], const int z_ind=-1) const;
	double dot3(const double aup[3], const double bup[3], const double xup[4], const int z_ind=-1) const;


public:

	virtual ~Grid() {}

	Transport* sim;

	std::string grid_type;

	Axis nu_grid_axis;
	vector<ScalarMultiDArray<NDIMS+1> > BB, abs_opac, scat_opac; // one for each species

	vector<SpectrumArray*> distribution;  // radiation energy density for each species in lab frame (erg/ccm. Integrated over bin frequency and direction)

	ScalarMultiDArray<NDIMS> alpha;
	ScalarMultiDArray<NDIMS> Q_annihil; // annihilation energy deposition rate (erg/ccm/s) (lab frame)
	ScalarMultiDArray<NDIMS> rho;       // density (g/cm^3)
	ScalarMultiDArray<NDIMS> T;         // gas temperature (K)
	ScalarMultiDArray<NDIMS> Ye;        // electron fraction
	ScalarMultiDArray<NDIMS> H_vis;     // specific heating rate (erg/s/g)

	MultiDArray<4,NDIMS> fourforce_abs, fourforce_emit;
	ScalarMultiDArray<NDIMS> l_abs, l_emit; // lepton number emission rate (cm^-3 s^-1) (comoving frame)

	// set everything up
	virtual void init(Lua* lua, Transport* insim);

	// write out zone information
	void         write_zones(const int iw) const;
	void         write_hdf5_data       (H5::H5File file) const;
	virtual void write_hdf5_coordinates(H5::H5File file) const=0;

	// radius using naieve coord transformation
	static double radius(const double xup[4]);
	virtual void axis_vector(vector<Axis>& axes) const = 0;

	// get directional indices from the zone index
	virtual void    zone_directional_indices(const int z_ind, vector<unsigned>& dir_ind) const=0;
	virtual void    dims                    (hsize_t dims[], const int size) const=0;
	virtual hsize_t dimensionality          () const=0;

	// describe zone
	virtual int    zone_index      (const double xup[4])                    const=0;
	virtual double zone_lab_3volume(const int z_ind)                        const=0;
	virtual double zone_lapse      (const int z_ind)                        const=0;
	virtual double zone_min_length (const int z_ind)                        const=0;
	virtual double zone_radius     (const int z_ind)                        const=0;
	virtual double zone_cell_dist  (const double p_xup[3], const int z_ind) const;
	virtual double zone_lorentz_factor(const int z_ind                    ) const=0;
	double         zone_com_3volume(const int z_ind)                        const;
	double         zone_4volume    (const int z_ind)                        const;
	double         zone_rest_mass  (const int z_ind)                        const;

	// global functions
	double total_rest_mass() const;

	// give the velocity vector at this point (PARTICLE COORDINATES)
	virtual void interpolate_fluid_velocity(const double xup[4], double v[3], int z_ind=-1) const = 0;
	void interpolate_metric(const double xup[4], Metric* g, const int z_ind=-1);

	// boundary conditions
	virtual void symmetry_boundaries(EinsteinHelper *eh, const double step_size) const=0;

	// move the particle
	void integrate_geodesic(EinsteinHelper *eh) const;

	// help with spawning particles
	void random_core_x_D(const double r_core, ThreadRNG *rangen, double xup[4], double D[3]) const;
	void isotropic_kup_tet(const double nu, double kup_tet[4], const double xup[4], ThreadRNG *rangen) const;
	void isotropic_direction(double D[3], ThreadRNG *rangen) const;
	virtual void sample_in_zone(const int z_ind, ThreadRNG *rangen, double xup[3]) const = 0;
};


#endif

