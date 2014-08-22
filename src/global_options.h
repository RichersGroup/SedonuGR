#ifndef _GLOBAL_OPTIONS_H
#define _GLOBAL_OPTIONS_H 1

#ifdef __INTEL_COMPILER
#pragma warning disable 161
#endif

#include <limits>
#include <cassert>
#include <mpi.h>
using real = float; // or float
const MPI_Datatype MPI_real = ( sizeof(real)==4 ? MPI_FLOAT : MPI_DOUBLE );
#define NaN std::numeric_limits<real>::quiet_NaN()
#define MAX std::numeric_limits<int>::max()

using namespace std;

#include "physical_constants.h"
namespace pc = physical_constants;

#endif

// NOTE: lots of things require double precision in random places because we're using
// physical units rather than non-dimensionalized units. It should be possible to turn
// some of the big data structures smaller, but with care
