#ifndef STATICVECTOR_HPP
#define STATICVECTOR_HPP

#include <iostream>
#include <cassert>
#include <cstdint>

namespace phys {

//Forward declaration of the staticVec class
template<class T, uint32_t N> class staticVec;

//Forward declaration of the output operator for staticVecs
template<class T, uint32_t N>
std::ostream& operator<<(std::ostream &os, const staticVec<T, N> &vec);

template<class T, uint32_t N>
class staticVec {
public:
	/* Default Constructor */
	staticVec(const T &a = 0);

	/* Special constructor for N == 2 */
	staticVec(const T &a, const T &b);

	/* Special constructor for N == 3 */
	staticVec(const T &a, const T &b, const T &c);

	/* Copy Constructor */
	staticVec(const staticVec &v);

	/* Copy Assignment Operator */
	const staticVec& operator=(const staticVec&);

	/* Assignment to scalar value */
	const staticVec& operator=(const T&);

	/* Read access */
	const T& operator[](uint32_t i) const;
	/* Write access */
	T& operator[](uint32_t i);

	/* Other member operators */
	staticVec& operator*=(const T&);
	staticVec& operator/=(const T&);

	staticVec& operator+=(const staticVec&);
	staticVec& operator-=(const staticVec&);

	/* Friend Functions */
	friend staticVec operator+(staticVec lhs, const staticVec &rhs) {
		lhs += rhs;
		return lhs;
	}

	friend staticVec operator-(staticVec lhs, const staticVec &rhs) {
		lhs -= rhs;
		return lhs;
	}

	friend staticVec operator*(staticVec lhs, const T &a) {
		lhs *= a;
		return lhs;
	}

	friend staticVec operator*(const T &a, staticVec rhs) {
		return rhs * a;
	}

	friend staticVec operator/(staticVec lhs, const T &a) {
		lhs /= a;
		return lhs;
	}

	friend std::ostream& operator<<(std::ostream &os,
			const staticVec<T, N> &v) {
		for (uint32_t i = 0; i < N; ++i) {
			os << v[i] << " ";
		}
		return os;
	}

private:
	T m_element[N];
};

using point2 = staticVec<double, 2>; //2D spatial coordinates
using point3 = staticVec<double, 3>; //3D spatial coordinates
using grid2 = staticVec<int, 2>; //2D discrete points
using grid3 = staticVec<int, 3>; //3D discrete points
using pixel2 = staticVec<uint32_t, 2>; //2D image points
using pixel3 = staticVec<uint32_t, 3>; //3D image points


} //namespace

#include "StaticVector.inl"

#endif //header guard
