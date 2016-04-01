#ifndef STATICVECTOR_HPP
#define STATICVECTOR_HPP

#include <iostream>
#include <cassert>

namespace phys{

	//Forward declaration of the staticVec class
	template<class T, size_t N> class staticVec;

	//Forward declaration of the output operator for staticVecs 
	template<class T, size_t N>
	std::ostream& operator<<(std::ostream& os, const staticVec<T,N>& vec); 

	template<class T, size_t N>
	class staticVec{
	public:
		/* Default Constructor */
		staticVec(	const T& a = 0	);

		/* Special constructor for N == 2 */
		staticVec(	const T& a, 
					const T& b	);

		/* Special constructor for N == 3 */
		staticVec(	const T& a, 
					const T& b,
					const T& c	);

		/* Copy Constructor */
		staticVec(	const staticVec& v );

		/* Copy Assignment Operator */
		const staticVec& operator=(const staticVec&);

		/* Assignment to scalar value */
		const staticVec& operator=(const T&);

		/* Read access */
		const T& operator[](size_t i) const;
		/* Write access */
		T& operator[](size_t i);

		/* Other member operators */
		staticVec& operator*=(const T&);
		staticVec& operator/=(const T&); 
		
		staticVec& operator+=(const staticVec&);
		staticVec& operator-=(const staticVec&); 

		/* Friend Functions */
		friend staticVec operator+(staticVec lhs, const staticVec& rhs)
		{
			lhs += rhs;
			return lhs;
		}

		friend staticVec operator-(staticVec lhs, const staticVec& rhs)
		{
			lhs -= rhs;
			return lhs;
		}

		friend staticVec operator*(staticVec lhs, const T& a)
		{
			lhs *= a;
			return lhs;
		}

		friend staticVec operator*(const T& a, staticVec rhs)
		{
			return rhs * a;
		}

		friend staticVec operator/(staticVec lhs, const T& a)
		{
			lhs /= a;
			return lhs;		
		}

		friend std::ostream& operator<< (std::ostream& os, const staticVec<T,N>& v)
		{
			for (size_t i = 0; i < N; ++i)
			{
				os << v[i] << " ";
			}
			return os;
		}

	private:
		T m_element[N];
	};

	typedef staticVec<double, 2> point2; //!< location within a continuous 2D space; could be cartesian (x,y), polar(r, theta), etc.
	typedef staticVec<double, 3> point3; //!< location witnin a continuous 3D space; (x,y,z), (r,phi,theta), etc.
	typedef staticVec<int, 2>	 grid2;  //!< location within a discrete 2D space; single plane atomic lattice for example.
	typedef staticVec<int, 3>	 grid3;  //!< location within a discrete 3D space; bulk atomic crystal for example.
	typedef staticVec<size_t, 2> pixel2; //!< location within a postive discrete 2D space; a greyscale digital image for example.
	typedef staticVec<size_t, 3> pixel3; //!< location within a postive discrete 3D space; an image with multiple channels for example.
}

#include "StaticVector.inl"

#endif
