#ifndef _SIMDARRAY_H
#define _SIMDARRAY_H 1

#include "global_options.h"
#include <iostream>
#include <array>
#include <cassert>
#include <limits>
#include <cstdalign>

using namespace std;


//===========//
// SIMDArray //
//===========//
template<typename T, unsigned len>
class SIMDArray : public array<T,len>{
public:

	//==============//
	// CONSTRUCTORS //
	//==============//
	SIMDArray(){}
	// construct from a SIMDArray
	template<typename Tin>
	explicit SIMDArray(const SIMDArray<Tin,len>& input){
		//#pragma omp simd
		for(unsigned i=0; i<len; i++)
			this->operator[](i) = input[i];
	}
	template<unsigned len2>
	SIMDArray(const SIMDArray<T,len2>& input){
		//#pragma omp simd
		for(unsigned i=0; i<len; i++)
			this->operator[](i) = i<len2 ? input[i] : NaN;
	}
	SIMDArray(const T& input){
		for(unsigned i=0; i<len; i++)
			this->operator[](i) = input;
	}
	// conditional construct
	SIMDArray(const SIMDArray<bool,len>& condition, const SIMDArray<T,len>& trueval, const SIMDArray<T,len>& falseval){
		if(      all(condition)) *this = trueval;
		else if(!any(condition)) *this = falseval;
		else{
			//#pragma omp simd
			for(unsigned i=0; i<len; i++){
				this->operator[](i) = condition[i] ? trueval[i] : falseval[i];
			}
		}
	}
	SIMDArray(const bool condition, const SIMDArray<T,len>& trueval, const SIMDArray<T,len>& falseval){
		*this = SIMDArray(SIMDArray<bool,len>(condition), trueval, falseval);
	}
	SIMDArray(const bool condition, const SIMDArray<T,len>& trueval, const T& falseval){
		*this = SIMDArray(SIMDArray<bool,len>(condition), trueval, SIMDArray<T,len>(falseval));
	}
	SIMDArray(const bool condition, const T& trueval, const SIMDArray<T,len>& falseval){
		*this = SIMDArray(SIMDArray<bool,len>(condition), SIMDArray<T,len>(trueval), falseval);
	}
	SIMDArray(const bool condition, const T& trueval, const T& falseval){
		*this = SIMDArray(SIMDArray<bool,len>(condition), SIMDArray<T,len>(trueval), SIMDArray<T,len>(falseval));
	}
	SIMDArray(const SIMDArray<bool,len>& condition, const SIMDArray<T,len>& trueval, const T& falseval){
		*this = SIMDArray(condition, trueval, SIMDArray<T,len>(falseval));
	}
	SIMDArray(const SIMDArray<bool,len>& condition, const T& trueval, const SIMDArray<T,len>& falseval){
		*this = SIMDArray(condition, SIMDArray<T,len>(trueval), falseval);
	}
	SIMDArray(const SIMDArray<bool,len>& condition, const T& trueval, const T& falseval){
		*this = SIMDArray(condition, SIMDArray<T,len>(trueval), SIMDArray<T,len>(falseval));
	}

	// ARITHMETIC OPERATORS
	 SIMDArray<T,len> operator%(const SIMDArray<T,len>& input) const{
		SIMDArray<T,len> result;
		//#pragma omp simd
		for(unsigned i=0; i<len; i++) result[i] = this->operator[](i) % input[i];
		return result;
	}
	SIMDArray<T,len> operator%(const T& input) const{
		return operator%(SIMDArray<T,len>(input));
	}
	SIMDArray<T,len> operator*(const SIMDArray<T,len>& input) const{
		SIMDArray<T,len> result;
		//#pragma omp simd
		for(unsigned i=0; i<len; i++) result[i] = this->operator[](i) * input[i];
		return result;
	}
	SIMDArray<T,len> operator*(const T& input) const{
		return operator*(SIMDArray<T,len>(input));
	}
	SIMDArray<T,len> operator+(const SIMDArray<T,len>& input) const{
		SIMDArray<T,len> result;;
		//#pragma omp simd
		for(unsigned i=0; i<len; i++) result[i] = this->operator[](i) + input[i];
		return result;
	}
	SIMDArray<T,len> operator+(const T& input) const{
		return operator+(SIMDArray<T,len>(input));
	}
	SIMDArray<T,len> operator/(const SIMDArray<T,len>& input) const{
		SIMDArray<T,len> result;
		//#pragma omp simd
		for(unsigned i=0; i<len; i++) result[i] = this->operator[](i) / input[i];
		return result;
	}
	SIMDArray<T,len> operator/(const T& input) const{
		return operator/(SIMDArray<T,len>(input));
	}
	SIMDArray<T,len> operator-(const SIMDArray<T,len>& input) const{
		SIMDArray<T,len> result;
		//#pragma omp simd
		for(unsigned i=0; i<len; i++) result[i] = this->operator[](i) - input[i];
		return result;
	}
	SIMDArray<T,len> operator-(const T& input) const{
		return operator-(SIMDArray<T,len>(input));
	}
	SIMDArray<T,len> operator-() const{
	  return (*this)*(-1);
	}

	// INCREMENT OPERATORS
	void operator+=(const SIMDArray<T,len>& input){
		*this = operator+(input);
	}
	void operator+=(const T& input){
		*this = operator+(input);
	}
	void operator-=(const SIMDArray<T,len>& input){
		*this = operator-(input);
	}
	void operator-=(const T& input){
		*this = operator-(input);
	}
	void operator*=(const SIMDArray<T,len>& input){
		*this = operator*(input);
	}
	void operator*=(const T& input){
		*this = operator*(input);
	}
	void operator/=(const SIMDArray<T,len>& input){
		*this = operator/(input);
	}
	void operator/=(const T& input){
		*this = operator/(input);
	}

	// COMPARISON OPERATORS
	SIMDArray<bool,len> operator==(const SIMDArray<T,len>& input) const{
		SIMDArray<bool,len> result = false;
		//#pragma omp simd
		for(unsigned i=0; i<len; i++) result[i] = (this->operator[](i) == input[i]);
		return result;
	}
	SIMDArray<bool,len> operator==(const T& input) const{
		return operator==(SIMDArray<T,len>(input));
	}
	SIMDArray<bool,len> operator!=(const SIMDArray<T,len>& input) const{
		SIMDArray<bool,len> result = false;
		//#pragma omp simd
		for(unsigned i=0; i<len; i++) result[i] = (this->operator[](i) != input[i]);
		return result;
	}
	SIMDArray<bool,len> operator!=(const T& input) const{
		return operator!=(SIMDArray<T,len>(input));
	}
	SIMDArray<bool,len> operator<=(const SIMDArray<T,len>& input) const{
		SIMDArray<bool,len> result = false;
		//#pragma omp simd
		for(unsigned i=0; i<len; i++) result[i] = (this->operator[](i) <= input[i]);
		return result;
	}
	SIMDArray<bool,len> operator<=(const T& input) const{
		return operator<=(SIMDArray<T,len>(input));
	}
	SIMDArray<bool,len> operator>=(const SIMDArray<T,len>& input) const{
		SIMDArray<bool,len> result = false;
		//#pragma omp simd
		for(unsigned i=0; i<len; i++) result[i] = (this->operator[](i) >= input[i]);
		return result;
	}
	SIMDArray<bool,len> operator>=(const T& input) const{
		return operator>=(SIMDArray<T,len>(input));
	}
	SIMDArray<bool,len> operator<(const SIMDArray<T,len>& input) const{
		SIMDArray<bool,len> result = false;
		//#pragma omp simd
		for(unsigned i=0; i<len; i++) result[i] = (this->operator[](i) < input[i]);
		return result;
	}
	SIMDArray<bool,len> operator<(const T& input) const{
		return operator<(SIMDArray<T,len>(input));
	}
	SIMDArray<bool,len> operator>(const SIMDArray<T,len>& input) const{
		SIMDArray<bool,len> result = false;
		//#pragma omp simd
		for(unsigned i=0; i<len; i++) result[i] = (this->operator[](i) > input[i]);
		return result;
	}
	SIMDArray<bool,len> operator>(const T& input) const{
		return operator>(SIMDArray<T,len>(input));
	}
	SIMDArray<bool,len> operator||(const SIMDArray<T,len>& input) const{
		SIMDArray<bool,len> result = false;
		//#pragma omp simd
		for(unsigned i=0; i<len; i++) result[i] = (this->operator[](i) || input[i]);
		return result;
	}
	SIMDArray<bool,len> operator||(const T& input) const{
		return operator||(SIMDArray<T,len>(input));
	}
	SIMDArray<bool,len> operator&&(const SIMDArray<T,len>& input) const{
		SIMDArray<bool,len> result = false;
		//#pragma omp simd
		for(unsigned i=0; i<len; i++) result[i] = (this->operator[](i) && input[i]);
		return result;
	}
	SIMDArray<bool,len> operator&&(const T& input) const{
		return operator&&(SIMDArray<T,len>(input));
	}
	SIMDArray<bool,len> operator&(const SIMDArray<T,len>& input) const{
		SIMDArray<bool,len> result = false;
		//#pragma omp simd
		for(unsigned i=0; i<len; i++) result[i] = (this->operator[](i) & input[i]);
		return result;
	}
	SIMDArray<bool,len> operator&(const T& input) const{
		return operator&(SIMDArray<T,len>(input));
	}
	SIMDArray<bool,len> operator<<(const SIMDArray<T,len>& input) const{
		SIMDArray<bool,len> result = false;
		//#pragma omp simd
		for(unsigned i=0; i<len; i++) result[i] = (this->operator[](i) << input[i]);
		return result;
	}
	SIMDArray<bool,len> operator<<(const T& input) const{
		return operator<<(SIMDArray<T,len>(input));
	}
	SIMDArray<bool,len> operator>>(const SIMDArray<T,len>& input) const{
		SIMDArray<bool,len> result = false;
		//#pragma omp simd
		for(unsigned i=0; i<len; i++) result[i] = (this->operator[](i) >> input[i]);
		return result;
	}
	SIMDArray<bool,len> operator>>(const T& input) const{
		return operator>>(SIMDArray<T,len>(input));
	}

	// masked modify
	void masked_set(const SIMDArray<bool,len>& mask, const SIMDArray<T,len>& trueval){
		if(all(mask)) *this = trueval;
		else{
			//#pragma omp simd
			for(unsigned i=0; i<len; i++){
				if(mask[i]) this->operator[](i) = trueval[i];
			}
		}
	}

	T sum(const SIMDArray<bool,len> mask){
		T result = 0;
		for(unsigned i=0; i<len; i++)
			if(mask[i]) result += this->operator[](i);
		return result;
	}
	T sum(){
		return sum(true);
	}

	T min(const SIMDArray<bool,len> mask){
		T result = this->operator[](0);
		for(unsigned i=1; i<len; i++) if(mask[i]){
			T next = this->operator[](i);
			result = next<result ? next : result;
		}
		return result;
	}
	T min(){
		return min(true);
	}
	T max(const SIMDArray<bool, len> mask){
		T result = this->operator[](0);
		for(unsigned i=1; i<len; i++) if(mask[i]){
			T next = this->operator[](i);
			result = next>result ? next : result;
		}
		return result;
	}
	T max(){
		return max(true);
	}

};

template<typename T, unsigned len>
SIMDArray<T,len> operator%(const T& input, const SIMDArray<T,len>& right){
	return SIMDArray<T,len>(input)%right;}
template<typename T, unsigned len>
SIMDArray<T,len> operator+(const T& input, const SIMDArray<T,len>& right){
	return SIMDArray<T,len>(input)+right;}
template<typename T, unsigned len>
SIMDArray<T,len> operator-(const T& input, const SIMDArray<T,len>& right){
	return SIMDArray<T,len>(input)-right;}
template<typename T, unsigned len>
SIMDArray<T,len> operator*(const T& input, const SIMDArray<T,len>& right){
	return SIMDArray<T,len>(input)*right;}
template<typename T, unsigned len>
SIMDArray<T,len> operator/(const T& input, const SIMDArray<T,len>& right){
	return SIMDArray<T,len>(input)/right;}
template<typename T, unsigned len>
SIMDArray<bool,len> operator==(const T& input, const SIMDArray<T,len>& right){
	return SIMDArray<T,len>(input)==right;}
template<typename T, unsigned len>
SIMDArray<bool,len> operator!=(const T& input, const SIMDArray<T,len>& right){
	return SIMDArray<T,len>(input)!=right;}
template<typename T, unsigned len>
SIMDArray<bool,len> operator>=(const T& input, const SIMDArray<T,len>& right){
	return SIMDArray<T,len>(input)>=right;}
template<typename T, unsigned len>
SIMDArray<bool,len> operator<=(const T& input, const SIMDArray<T,len>& right){
	return SIMDArray<T,len>(input)<=right;}
template<typename T, unsigned len>
SIMDArray<bool,len> operator>(const T& input, const SIMDArray<T,len>& right){
	return SIMDArray<T,len>(input)>right;}
template<typename T, unsigned len>
SIMDArray<bool,len> operator<(const T& input, const SIMDArray<T,len>& right){
	return SIMDArray<T,len>(input)<right;}
template<typename T, unsigned len>
SIMDArray<bool,len> operator&&(const T& input, const SIMDArray<T,len>& right){
	return SIMDArray<T,len>(input)&&right;}
template<typename T, unsigned len>
SIMDArray<bool,len> operator||(const T& input, const SIMDArray<T,len>& right){
	return SIMDArray<T,len>(input)||right;}
template<typename T, unsigned len>
SIMDArray<bool,len> operator&(const T& input, const SIMDArray<T,len>& right){
	return SIMDArray<T,len>(input)&right;}
template<typename T, unsigned len>
SIMDArray<bool,len> operator<<(const T& input, const SIMDArray<T,len>& right){
	return SIMDArray<T,len>(input)<<right;}
template<typename T, unsigned len>
SIMDArray<bool,len> operator>>(const T& input, const SIMDArray<T,len>& right){
	return SIMDArray<T,len>(input)>>right;}

// ANY / ALL
template<unsigned len>
bool any(const SIMDArray<bool,len>& input){
	bool result = false;
	for(unsigned i=0; i<len; i++) result = (result || input[i]);
	return result;
}
template<unsigned len>
bool all(const SIMDArray<bool,len>& input){
	bool result = true;
	for(unsigned i=0; i<len; i++) result = (result && input[i]);
	return result;
}

// MISCELLANEOUS
template<typename T, unsigned n>
std::ostream& operator<<(std::ostream& out, const SIMDArray<T,n>& SIMDArray){
	out << "[ ";
	for(unsigned i=0; i<n; i++) out << SIMDArray[i] << " ";
	out << "]";
    return out;
}

template<unsigned len>
SIMDArray<bool,len> operator!(const SIMDArray<bool,len>& ch){
	SIMDArray<bool,len> result;
	//#pragma omp simd
	for(unsigned i=0; i<len; i++)
		result[i] = !ch[i];
	return result;
}
template<typename T, unsigned len>
SIMDArray<T,len> sqrt(const SIMDArray<T,len>& ch){
	SIMDArray<T,len> result;
	//#pragma omp simd
	for(unsigned i=0; i<len; i++)
		result[i] = sqrt(ch[i]);
	return result;
}

template<typename T, unsigned len>
SIMDArray<T,len> abs(const SIMDArray<T,len>& ch){
	SIMDArray<T,len> result;
	//#pragma omp simd
	for(unsigned i=0; i<len; i++)
		result[i] = abs(ch[i]);
	return result;
}
template<typename T, unsigned len>
SIMDArray<T,len> log(const SIMDArray<T,len>& ch){
	SIMDArray<T,len> result;
	//#pragma omp simd
	for(unsigned i=0; i<len; i++)
		result[i] = log(ch[i]);
	return result;
}
template<typename T, unsigned len>
SIMDArray<T,len> exp(const SIMDArray<T,len>& ch){
	SIMDArray<T,len> result;
	//#pragma omp simd
	for(unsigned i=0; i<len; i++)
		result[i] = exp(ch[i]);
	return result;
}
template<typename T, unsigned len, typename T2>
SIMDArray<T,len> pow(const SIMDArray<T,len>& ch, const T2 exponent){
	SIMDArray<T,len> result;
	//#pragma omp simd
	for(unsigned i=0; i<len; i++)
		result[i] = pow(ch[i], exponent);
	return result;
}
template<typename T, unsigned len>
SIMDArray<T,len> min(const SIMDArray<T,len>& in1, const SIMDArray<T,len>& in2){
	SIMDArray<T,len> result;
	//#pragma omp simd
	for(unsigned i=0; i<len; i++)
		result[i] = min(in1[i], in2[i]);
	return result;
}
template<typename T, unsigned len>
SIMDArray<T,len> min(const SIMDArray<T,len>& in1, const T& in2){
	return min(in1, SIMDArray<T,len>(in2));
}
template<typename T, unsigned len>
SIMDArray<T,len> min(const T& in2, const SIMDArray<T,len>& in1){
	return min(in1, SIMDArray<T,len>(in2));
}
template<typename T, unsigned len>
SIMDArray<T,len> max(const SIMDArray<T,len>& in1, const SIMDArray<T,len>& in2){
	SIMDArray<T,len> result;
	//#pragma omp simd
	for(unsigned i=0; i<len; i++)
		result[i] = max(in1[i], in2[i]);
	return result;
}
template<typename T, unsigned len>
SIMDArray<T,len> max(const SIMDArray<T,len>& in1, const T& in2){
	return max(in1, SIMDArray<T,len>(in2));
}
template<typename T, unsigned len>
SIMDArray<T,len> max(const T& in2, const SIMDArray<T,len>& in1){
	return max(in1, SIMDArray<T,len>(in2));
}
template<typename T, unsigned len>
SIMDArray<T,len> atan2(const SIMDArray<T,len>& in1, const SIMDArray<T,len>& in2){
	SIMDArray<T,len> result;
	//#pragma omp simd
	for(unsigned i=0; i<len; i++)
		result[i] = atan2(in1[i], in2[i]);
	return result;
}
template<typename T, unsigned len>
SIMDArray<T,len> cos(const SIMDArray<T,len>& in1){
	SIMDArray<T,len> result;
	//#pragma omp simd
	for(unsigned i=0; i<len; i++)
		result[i] = cos(in1[i]);
	return result;
}
template<typename T, unsigned len>
SIMDArray<T,len> sin(const SIMDArray<T,len>& in1){
	SIMDArray<T,len> result;
	//#pragma omp simd
	for(unsigned i=0; i<len; i++)
		result[i] = sin(in1[i]);
	return result;
}

template<unsigned len>
bool count(const SIMDArray<bool,len> input){
	unsigned result = 0;
	for(unsigned i=0; i<len; i++) result += (unsigned)input[i];
	return result;
}

template<typename T, unsigned len>
double radius(const SIMDArray<T,len>& x){
	assert(len >= 3);
	return sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
}
template<unsigned len, unsigned chunksize>
SIMDArray<double,chunksize> radius(const SIMDArray<SIMDArray<double,chunksize>, len>& x){
	assert(len >= 3);
	return sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
}

#endif
