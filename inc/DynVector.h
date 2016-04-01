#ifndef PHYSVECTOR_HPP
#define PHYSVECTOR_HPP

#include "Complex.h"
#include "Helpers.h"
#include <cassert>


namespace phys{
	/*****************************************************************
	*	Overload operators for the std::vector class inside our namespace
		@TODO: add global operators to call these so that we don't have
		to use using namespace, or the function-type call
	*****************************************************************/

	/* Increment the elements in the lhs by the corresponding elements in the rhs
		- returns a non-const reference for chaining, e.g. w = u += v; */
	template<typename T>
	std::vector<T>& operator+=(std::vector<T>& lhs, const std::vector<T>&rhs);

	/* Decrement the elements in the lhs by the corresponding elements in the rhs
	- returns a non-const reference for chaining, e.g. w = u -= v; */
	template<typename T>
	std::vector<T>& operator-=(std::vector<T>& lhs, const std::vector<T>&rhs);

	template<typename T>
	std::vector<T>& operator+=(std::vector<T>& lhs, const T& rhs);

	template<typename T>
	std::vector<T>& operator-=(std::vector<T>& lhs, const T& rhs);

	/* Multiplies the vector by a scalar */
	template<typename T>
	std::vector<T>& operator*=(std::vector<T>& lhs, const T& rhs);

	/* Divides the vector by a scalar */
	template<typename T>
	std::vector<T>& operator/=(std::vector<T>& lhs, const T& rhs);

	//addition of vectors - corresponding elements are added
	template<typename T>
	const std::vector<T> operator+(const std::vector<T>& lhs, const std::vector<T>& rhs);

	//subtraction of vector - corresponding elements are subtracted
	template<typename T>
	const std::vector<T> operator-(const std::vector<T>& lhs, const std::vector<T>& rhs);

	//Vector inner product
	template<typename T>
	const T operator*(const std::vector<T>& lhs, const std::vector<T>& rhs);

	//squared norm of vector u.
	template<typename T>
	const T squaredNorm(const std::vector<T>& v);

	//compute the norm of a vector, default value for p is 2 == Euclidean norm
	template<typename T>
	const T norm(const std::vector<T>&v, int p = 2);

	//postive of a vector
	template<typename T>
	const std::vector<T> operator+(const std::vector<T>& v);

	//negate a vector
	template<typename T>
	const std::vector<T> operator-(const std::vector<T>& v);

	//vector times scalar
	template<typename T>
	const std::vector<T> operator*(const std::vector<T>& lhs, const T& rhs);

	//scalar times vector
	template<typename T>
	const std::vector<T> operator*(const T& lhs, const std::vector<T>& rhs );

	//vector divide by scalar
	template<typename T>
	const std::vector<T> operator/(const std::vector<T>& lhs, const T& rhs);

	//Output operator for vectors - returns non-const reference for chaining
	template<typename T>
	std::ostream & operator<<(std::ostream& os, const std::vector<T>& v);

	/* Read in a vector from a text file: use std::getline to get string argument.*/
	template<typename T>
	const std::vector<T> read_vector(const std::string& line, const T& t = T());

	/* select a contiguous sub-vector in the range [start,end] (inclusive of end) */
	template<typename T>
	const std::vector<T> sub_vector( const std::vector<T>& vec, size_t start, size_t end); 

	/* Type definitions for convenience */
	typedef std::vector<double>		stdVec_d;
	typedef std::vector<int>		stdVec_i;
	typedef std::vector<size_t>		stdVec_s;

	typedef std::vector<phys::complex>	stdVec_c;
}

template<typename T>
std::ostream& operator<< (std::ostream& os, const std::vector<T>& v)
{
	return phys::operator<<(os, v);
}//global output stream operator for covenience

#include "DynVector.inl"

#endif
