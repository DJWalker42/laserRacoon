#ifndef PHYSCOMPLEX_HPP
#define PHYSCOMPLEX_HPP

#include "StaticVector.h"

namespace phys{

	class complex : public point2{
	public:
		complex();
		complex(const point2& p);
		complex(const complex& c);
		complex(double re, double im);

		//explicitly define assignment to double
		const complex& operator=(double);

		//We redefine these operators as they have different meaning from the base class
		const complex& operator+= (double);
		const complex& operator-= (double);
		
		//this operator is not defined for the base class but has meaning for the complex class.
		const complex& operator*= (const complex&);
		const complex& operator/= (const complex&);

		const complex operator/=(double); 

		//complex conjugate using the unary + operator
		friend complex operator+(const complex& c);

		void conjugate();

		double norm() const;
		double argument() const;
		complex power(double n) const; 
	};

	const complex operator+ (double, const complex&);
	const complex operator+ (const complex&, double);
	const complex operator- (double, const complex&);
	const complex operator- (const complex&, double);
	const complex operator* (const complex&, const complex&);
	const complex operator/ (const complex&, const complex&);

	double norm( const complex& );
	complex power( const complex&, double n);
}
#endif
