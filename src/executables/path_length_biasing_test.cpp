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
#include <algorithm>
#include "Transport.h"
#include "global_options.h"

using namespace std;

int main(int argc, char **argv)
{
	// set up MPI
	MPI_Init( &argc, &argv );

	// set up the transport module (includes the grid)
	int nbins = 0, nsamples = 0,niter = 0;
	double max_tau_bin = 0, abs_frac=0;
	Transport sim;

	// setup and seed random number generator(s)
	sim.rangen.init();

	// read in parameters from the command line
	sscanf(argv[1], "%lf",  &abs_frac);
	sscanf(argv[2], "%d",  &nbins);
	sscanf(argv[3], "%lf", &max_tau_bin);
	sscanf(argv[4], "%d",  &nsamples);
	sscanf(argv[5], "%d",  &niter);
	sim.bias_path_length=1;
	sim.max_path_length_boost = INFINITY;
	sim.min_packet_energy = 1e-6;

	// set up the data vectors
	vector<double> grid             = vector<double>(nbins,0);
	vector<double> energy           = vector<double>(nbins,0);
	vector<double> expected_energy  = vector<double>(nbins,0);
	vector<int>    packets          = vector<int>(nbins,0);
	vector<double> expected_packets = vector<double>(nbins,0);
	double alpha = 1./(1.0-abs_frac);
	for(int i=0; i<(int)grid.size(); i++){
		grid[i] = (i+1)*max_tau_bin / (double)nbins;
		expected_energy[i]  = ( i==0 ? 1 : exp(-grid[i-1])       ) - exp(-grid[i]      );
		expected_packets[i] = ( i==0 ? 1 : exp(-grid[i-1]/alpha) ) - exp(-grid[i]/alpha);
	}

	// set up the particle
	double v[3] = {0,0,0};
	LorentzHelper lh(false);
	lh.set_v(v,3);
	int n_above = 0;
	double e_above = 0;

	double avg_e = 0;
	int n_zero = 0;
	int n_dead = 0;

	// sample lots of particles
	for(int i=0; i<nsamples; i++){
		{
			Particle p;
			p.e = 1.0;
			p.fate = moving;
			double opacity = 1.0;//sim.rangen.uniform();
			lh.set_p<com>(&p);
			lh.set_opac<com>(opacity*abs_frac, opacity*(1.0-abs_frac));
		}
		for(int j=0; j<niter; j++) if(lh.p_fate()==moving){
			sim.sample_tau(&lh);
			while(lh.p_e(com)<=sim.min_packet_energy && lh.p_fate()==moving){
				if(sim.rangen.uniform() < 0.5) lh.set_p_fate(rouletted);
				else lh.scale_p_e(2.0);
			}
		}

		// log the sample
		int ind = upper_bound(grid.begin(), grid.end(), lh.p_tau()) - grid.begin();
		if(ind<(int)grid.size()){
			packets[ind]++;
			energy[ind] += lh.p_e(com);
			avg_e += lh.p_e(com);
		}
		else{
			n_above ++;
			e_above += lh.p_e(com);
		}
		if(lh.p_e(com)==0) n_zero++;
		if(lh.p_fate()!=moving) n_dead++;
		else PRINT_ASSERT(lh.p_e(com), >=, sim.min_packet_energy);
	}

	// calculate energy variance
	double delta2 = 0;
	double weight_sum = 0;
	for(int i=0; i<(int)grid.size(); i++){
		double tmp = (energy[i]/(double)nsamples-expected_energy[i]);
		delta2 += expected_energy[i] * tmp*tmp;
		weight_sum += expected_energy[i];
	}
	delta2 /= weight_sum;

	// print the results
	cout << "variance = " << sqrt(delta2) << endl;
	cout << "n_above = " << (double)n_above/(double)nsamples << "%" << endl;
	cout << "e_above = " << e_above/(double)nsamples << "%" << endl;
	cout << "average e = " << avg_e/(double)nsamples << endl;
	cout << "n_zero = " << (double)n_zero/(double)nsamples << endl;
	cout << "n_dead = " << (double)n_dead/(double)nsamples << endl;
	cout << "1-tau\t2-Ebar\t3-Nbar\t4-expected_energy\t5-expected_packets" << endl;
	for(int i=0; i<(int)grid.size(); i++){
		cout << grid[i] << '\t';
		cout << energy[i]          / (double)nsamples << '\t';
		cout << (double)packets[i] / (double)nsamples << '\t';
		cout << expected_energy[i] << '\t';
		cout << expected_packets[i] << endl;
	}

}

