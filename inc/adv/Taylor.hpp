#include <physVector.hpp>
#include <adv/polynomial.hpp>
#include <physMaths.hpp>

namespace phys{

	template<class T>
	const std::vector<T> TaylorScheme(	const std::vector<T>& u0,
									const matrix<T>& S,
									const list<std::vector<T>>& f, 
									const T& h )
	{
		T h_i_fact = 1;
		std::vector<T> sum(u0.size());
		std::vector<T> uDerivative = u0;
		for(size_t i = 0; i < f.size(); i++)
		{
			sum += h_i_fact  * uDerivative;
			uDerivative = S * uDerivative + f[i];
			h_i_fact  *= h/(i+1);
		}
		return sum;
	}

	template<class T>
	const std::vector<T> error(	const std::vector<T>&boundForU,
							const matrix<T>&S,
							const list<std::vector<T> >&f,
							const T& h)
	{
		T h_dash;
		h_dash = 1;
		std::vector<T> uDerivative = boundForU;
		for(size_t i = 0; i < f.size(); i++)
		{
			uDerivative = S * uDerivative + f[i];
			h_dash *= h/(i+1);
		}
		return h_dash * uDerivative;
	} // error in Taylor scheme

	template<class T>
	void deriveKS(	const T&c, 
					const T&r,
					list<T>&u, 
					list<T>&v,
					list<T>&w)
	{
		list<T> Rinverse = deriveRinverse(r, u.size());
		for(size_t i = 0; i < u.size()-1; i++)
		{
			u(i+1) = v[i]-deriveProduct(u,Rinverse,i);
			v(i+1) = w[i];
			w(i+1) = (-0.5)*deriveProduct(u,u,i)
				- deriveProduct(w,Rinverse,i)
				- v[i] + (i ? 0. : (c*c) );
		}
	} // derivatives in KS equation

	template<class T>
	const std::vector<T> TaylorKS(	const T&c, 
								const T&r, 
								const T&h, 
								size_t n,
								const std::vector<T>&u0,
								const std::vector<T>&bound)
	{
		list<T> u(n,0.);
		list<T> v(n,0.);
		list<T> w(n,0.);
		u(0) = u0[0];
		v(0) = u0[1];
		w(0) = u0[2];
		deriveKS(c,r,u,v,w);
		std::vector<T> result(3);
		result[0] = HornerTaylor(u,h);
		result[1] = HornerTaylor(v,h); 
		result[2] = HornerTaylor(w,h));
		u(0) = bound[0];
		v(0) = bound[1];
		w(0) = bound[2];
		deriveKS(c,r,u,v,w);
		std::vector<T> highDerivative(3);
		highDerivative[0] = u[n-1]; 
		highDerivative[1] = v[n-1]; 
		highDerivative[2] = w[n-1];
		std::vector<T> error = (power(h,n-1)/factorial(n-1)) * highDerivative;
		return result + error;
	} // Taylor for KS equation + error estimate


}