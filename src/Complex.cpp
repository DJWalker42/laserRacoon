#include <cmath>

#include "Complex.h"

namespace phys{

	/***************************************************
	* Constructors
	**************************************************/
	complex::complex() : point2() 
	{}

	complex::complex(const point2& p) : point2(p)
	{}

	complex::complex(const complex& c) : point2(c[0], c[1])
	{}

	complex::complex(double re, double im) : point2(re, im)
	{}

	/***************************************************
	* Member operators
	***************************************************/
	const complex& complex::operator= (double a)
	{
		(*this)[0] = a;
		(*this)[1] = a;
		return *this;
	}// assign to a double value.

	const complex& complex::operator+= (double a)
	{
		(*this)[0] += a;
		return *this;
	} // adding a real number

	const complex& complex::operator-= (double a)
	{
		(*this)[0] -= a;
		return *this;
	}// subtracting a real number

	const complex& complex::operator*=(const complex& c)
	{
		double keep = (*this)[0];
		(*this)[0] = keep * c[0] - (*this)[1] * c[1]; 
		(*this)[1] = keep * c[1] + (*this)[1] * c[0];
		return *this;
	}// complex number multiplication


	const complex& complex::operator/=(const complex& c)
	{
		return *this *= (+c) / squaredNorm(c); 
	} // complex number division

	const complex complex::operator/=(double d)
	{
		return *this / d;
	}

	/**********************************************************
	*	Member functions
	**********************************************************/
	void complex::conjugate()
	{
		(*this)[1] = -(*this)[1];
	}

	double complex::norm() const
	{
		return std::sqrt((*this)[0] * (*this)[0] + (*this)[1] * (*this)[1]);
	}

	double complex::argument() const
	{
		return atan2( (*this)[1], (*this)[0] );
	}

	complex complex::power(double n) const
	{
		double arg = n * argument();
		return complex(cos(arg), sin(arg)) * pow(norm(), n);
	}


	/**********************************************************
	*	Non-Member operators
	**********************************************************/
	const complex operator+ (double a, const complex& c)
	{
		return complex(c) += a;
	}// real plus complex

	const complex operator+ ( const complex& c, double a)
	{
		return complex(c) += a;
	}// complex plus real

	const complex operator- (double a, const complex& c)
	{
		return complex(a, 0.) - c;
	}// real minus complex

	const complex operator- (const complex& c, double a)
	{
		return complex(c) -= a;
	}// complex minus real

	const complex operator* (const complex& c, const complex& d)
	{
		return complex(c) *= d; 
	}

	const complex operator/ (const complex& c, const complex& d)
	{
		return complex(c) /= d;
	}

	double norm (const complex& c)
	{
		return sqrt(c[0] * c[0] + c[1] * c[1]); 
	}

	complex power(const complex&c, double n)
	{
		return c.power(n);
	}

	//friend complex conjugate operator+
	complex operator+(const complex& c)
	{
		return complex(c[0], -c[1]);
	}

}
