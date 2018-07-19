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
#include <cmath>
#include <atomic>
#include <array>

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


#if DEBUG==1
#define PRINT_ASSERT(a,op,b)                         \
do {                                                 \
	if(!((a) op (b))) std::cout << (a) << " " << (b) << std::endl; \
	assert(a op b);                                  \
} while (0)
#else
#define PRINT_ASSERT(a,op,b) do{} while(false)
#endif

template<typename T>
class ATOMIC : public std::atomic<T>
{
 public:
  ATOMIC() = default;

  constexpr ATOMIC(T desired) :
  std::atomic<T>(desired)
    {}

  constexpr ATOMIC(const ATOMIC<T>& other) :
  ATOMIC(other.load(std::memory_order_relaxed))
    {}

  ATOMIC& operator=(const ATOMIC<T>& other) {
    this->store(other.load(std::memory_order_relaxed), std::memory_order_relaxed);
    return *this;
    }

  void operator+=(T in) {
    T current = this->load();
    while (!this->compare_exchange_weak(current, current + in, std::memory_order_relaxed));
  }
  void operator-=(T in) {
    T current = this->load();
    while (!this->compare_exchange_weak(current, current - in, std::memory_order_relaxed));
  }
  template<typename Tin>
  void operator*=(Tin in) {
    operator=((*this)*in);
  }
};


//=======//
// TUPLE //
//=======//
template<typename T, size_t len>
class Tuple : public std::array<T,len>{
public:

	Tuple(){}

	Tuple(T in){
		for(size_t i=0; i<len; i++) this->operator[](i) = in;
	}
	template<typename Tin>
	  Tuple(Tuple<Tin,len> in){
		for(size_t i=0; i<len; i++) this->operator[](i) = in[i];
	}

  inline size_t size() const{return len;}

  template<typename Tin>
  Tuple<T,len>& operator=(const Tuple<Tin,len>& input){
	  for(size_t i=0; i<len; i++) this->operator[](i) = input[i];
	  return *this;
  }
  Tuple<T,len>& operator=(const double input){
	  for(size_t i=0; i<len; i++) this->operator[](i) = input;
	  return *this;
  }
  const Tuple<T,len> operator*(const double scale) const{
	  Tuple<T,len> result;// = *this;
	  for(size_t i=0; i<len; i++) result[i] = this->operator[](i) * scale;
	  return result;
  }
  const Tuple<T,len> operator/(const double scale) const{
	  double inv_scale = 1./scale;
	  return operator*(inv_scale);
  }
  template<typename Tin>
  const Tuple<T,len> operator+(const Tuple<Tin,len>& input) const{
	  Tuple<T,len> result;// = *this;
	  for(size_t i=0; i<len; i++) result[i] = this->operator[](i) + input[i];
	  return result;
  }
  template<typename Tin>
  const Tuple<T,len> operator-(const Tuple<Tin,len>& input) const{
	  Tuple<T,len> result;// = *this;
	  for(size_t i=0; i<len; i++) result[i] = this->operator[](i) - input[i];
	  return result;
  }
  template<typename Tin>
  void operator+=(const Tuple<Tin,len>& input){
	  *this = (*this) + input;
  }
  void operator*=(const double scale){
	  *this = (*this) * scale;
  }
  template<typename Tin>
  bool operator==(const Tuple<Tin,len>& input){
	  bool isequal = true;
	  for(size_t i=0; i<len; i++) isequal = isequal && (this->operator[](i) == input[i]);
	  return isequal;
  }
};

template<typename T, size_t n>
inline std::ostream& operator<<(std::ostream& out, const Tuple<T,n>& tuple){
	out << "[ ";
	for(size_t i=0; i<n; i++) out << tuple[n] << " ";
	out << "]";
    return out;
}

template<size_t n>
inline double radius(const Tuple<double,n>& x){
	PRINT_ASSERT(n,>=,3);
	return sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
}
inline double radius(const double x[3]){
	return sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
}

#endif

// NOTE: lots of things require double precision in random places because we're using
// physical units rather than non-dimensionalized units. It should be possible to turn
// some of the big data structures smaller, but with care
