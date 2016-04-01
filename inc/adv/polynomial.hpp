#ifndef POLYNOMIAL_HPP
#define POLYNOMIAL_HPP
#include <adv/physList.hpp>
#include <PascalTri.hpp>


namespace phys{

	template<class T> 
	class Polynomial : public list<T>{
	public:
		/*	construct an n-1 degree Polynomial (n terms) with zero coefficients */
		Polynomial(size_t n=0);

		/*	construct an n-1 degree Polynomial (n terms), each with the coefficient a */
		Polynomial(size_t n, const T&a);

		~Polynomial()
		{} // destructor -- deletion of pointers handled by ~list()

		size_t degree() const
		{
			return number-1;
		} // degree of Polynomial
	};

	// compound addtion of polynomials
	template<class T>
	const Polynomial<T>&operator+=(Polynomial<T>& p, const Polynomial<T>&q);

	// add two polynomials
	template<class T>
	const Polynomial<T>operator+(const Polynomial<T>& p,
			const Polynomial<T>&q);

	// Polynomial by scalar
	template<class T>
	const Polynomial<T>&operator*=(Polynomial<T>& p, const T&a);

	// scalar times Polynomial
	template<class T>
	const Polynomial<T>operator*(const T&a, const Polynomial<T>&p);

	// compound multiply by Polynomial
	template<class T>
	Polynomial<T>&operator*=(Polynomial<T>&p, const Polynomial<T>&q);

	// multiply two polynomials
	template<class T>
	Polynomial<T>operator*(const Polynomial<T>&p,const Polynomial<T>&q);

	// calculate a Polynomial at point x
	template<class T>
	const T calculatePolynomial(const Polynomial<T>&p, const T&x);

	// Horner algorithm to calculate a Polynomial
	template<class T>
	const T HornerPolynomial(const Polynomial<T>&p, const T&x);

	// derivatives of 1/r
	template<class T>
	const list<T>deriveRinverse(const T&r, size_t n);

	// Taylor approximation to f(x+h)
	template<class T>
	const T Taylor(const list<T>&f, const T&h);

	// Horner algorithm for Taylor approximation
	template<class T>
	const T HornerTaylor(const list<T>&f, const T&h);

	// nth derivative of a product
	template<class T>
	const T deriveProduct(const list<T>&f, const list<T>&g, size_t n);

	// integral in the unit interval
	template<class T>
	const T integral(const Polynomial<T>&p);

	// integral in the triangle
	template<class T>
	const T integral(const Polynomial<Polynomial<T> >&p);
}

#include <adv/Polynomial.inl>

#endif