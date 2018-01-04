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

#ifndef _GLOBAL_OPTIONS_H
#define _GLOBAL_OPTIONS_H 1

#ifdef __INTEL_COMPILER
#pragma warning disable 161
#endif

#include <limits>
#include <cassert>
#include <string>
#include <iostream>
#include "physical_constants.h"
#include "H5Cpp.h"

#define NDIMS 1

//using real = float; // or float
//const MPI_Datatype MPI_real = ( sizeof(real)==4 ? MPI_FLOAT : MPI_DOUBLE );
#define NaN std::numeric_limits<double>::quiet_NaN()
#define MAXLIM std::numeric_limits<int>::max()
#define TINY 1e-5

inline std::string trim(const std::string s)
{
	std::string trimmed = s;
	std::string::size_type pos = trimmed.find_last_not_of(' ');
	if(pos != std::string::npos)
	{
		if (trimmed.length()!=pos+1)//if there are trailing whitespaces erase them
			trimmed.erase(pos+1);
		pos = trimmed.find_first_not_of(' ');
		if(pos!=0) //if there are leading whitespaces erase them
			trimmed.erase(0, pos);
	}
	else trimmed="";
	return trimmed;
}

inline bool hdf5_dataset_exists(const char* filename, const char* datasetname){
	bool exists = true;

	// Temporarily turn off error printing
	H5E_auto2_t func;
	void* client_data;
	H5::Exception::getAutoPrint(func,&client_data);
	H5::Exception::dontPrint();

	// See if dataset exists
	H5::H5File file(filename, H5F_ACC_RDONLY);
	H5::DataSet dataset;
	try{
		dataset = file.openDataSet(datasetname);
	}
	catch(H5::FileIException& exception){
		exists = false;
	}

	// Turn error printing back on
	H5::Exception::setAutoPrint(func,client_data);
	file.close();

	return exists;
}


#ifndef NDEBUG
#define PRINT_ASSERT(a,op,b)                         \
do {                                                 \
	if(!(a op b)) std::cout << (a) << " " << (b) << std::endl; \
	assert(a op b);                                  \
} while (0)
#else
#define PRINT_ASSERT(a,op,b)
#endif

//=======//
// TUPLE //
//=======//
template<typename T, unsigned len>
class Tuple{
public:
  T vals[len];

  unsigned int size() const{return size;}

  const T& operator[](const unsigned int i) const {return vals[i];}
  T& operator[](const unsigned int i){return vals[i];}
  Tuple<T,len>& operator=(const Tuple<T,len> input){
	  for(unsigned i=0; i<len; i++) this->vals[i] = input.vals[i];
	  return *this;
  }
  Tuple<T,len>& operator=(const double input){
	  for(unsigned i=0; i<len; i++) this->vals[i] = input;
	  return *this;
  }
  const Tuple<T,len> operator*(const double scale) const{
	  Tuple<T,len> result = *this;
	  for(unsigned i=0; i<len; i++) result.vals[i] *= scale;
	  return result;
  }
  const Tuple<T,len> operator/(const double scale) const{
	  double inv_scale = 1./scale;
	  return operator*(inv_scale);
  }
  const Tuple<T,len> operator+(const Tuple<T,len> input) const{
	  Tuple<T,len> result = *this;
	  for(unsigned i=0; i<len; i++) result.vals[i] += input.vals[i];
	  return result;
  }
  const Tuple<T,len> operator-(const Tuple<T,len> input) const{
	  Tuple<T,len> result = *this;
	  for(unsigned i=0; i<len; i++) result.vals[i] -= input.vals[i];
	  return result;
  }
  Tuple<T,len>& operator+=(const Tuple<T,len> input){
	  for(unsigned i=0; i<len; i++) this->vals[i] += input.vals[i];
	  return *this;
  }
  Tuple<T,len>& operator*=(const double scale){
	  for(unsigned i=0; i<len; i++) this->vals[i] *= scale;
	  return *this;
  }
};


#endif

// NOTE: lots of things require double precision in random places because we're using
// physical units rather than non-dimensionalized units. It should be possible to turn
// some of the big data structures smaller, but with care
