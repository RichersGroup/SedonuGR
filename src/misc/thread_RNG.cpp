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

#include <ctime>
#include <mpi.h>
#include <omp.h>
#include "thread_RNG.h"

//-----------------------------------------------------------------
// initialize the RNG system
//-----------------------------------------------------------------
// ASSUMES the number of threads remains constant so it only has to be initialized once
void thread_RNG::init(){
	int my_mpiID;
	MPI_Comm_rank(MPI_COMM_WORLD, &my_mpiID);

	// set up the stuff that creates the random number generators
	const gsl_rng_type* TypeR = gsl_rng_default;
	gsl_rng_env_setup();

    #pragma omp parallel default(none) shared(my_mpiID,TypeR,gsl_rng_default_seed)
	{
        #ifdef _OPENMP
		const int nthreads = omp_get_num_threads();
		const int my_ompID = omp_get_thread_num();
		#else
		const int nthreads = 1;
		const int my_ompID = 0;
		#endif

		// assign a unique RNG to each thread
        #pragma omp single
		generators.resize(nthreads);
		for(int i=0; i<nthreads; i++){
			gsl_rng_default_seed = (unsigned int)time(NULL) + my_mpiID*nthreads + my_ompID;
			generators[i] = gsl_rng_alloc (TypeR);
		}
	}
}



//-----------------------------------------------------------------
// return a uniformily distributed random number (thread safe)
//-----------------------------------------------------------------
double thread_RNG::uniform(){
    #ifdef _OPENMP
	const int my_ompID = omp_get_thread_num();
    #else
	const int my_ompID = 0;
    #endif

	return gsl_rng_uniform(generators[my_ompID]);
}
