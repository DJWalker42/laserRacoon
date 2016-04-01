#include <cassert>

#include "DataFit.h"
#include "Maths.h"


namespace phys{
	namespace dfit{

		/** Linear least square fit to the data 
		*/
		std::pair<double,double> llsq( const stdVec_d& x, const stdVec_d& y )
		{
			assert(x.size() == y.size());

			//quick return if called with empty parameters
			if(x.empty() || y.empty()){
				return std::pair<double,double> (0.0, 0.0);	
			}		
			size_t n = x.size();
			std::pair<double,double> retval;
			//  Special cases.
			if ( n == 1 ){
				retval.first = 0.0;
				retval.second = y[0];
			} else if ( n == 2 ) {
				double m = (y[1] - y[0])/(x[1] - x[0]);
				double c = y[0] - m * x[0];
				retval.first = m;
				retval.second = c;
			} else {
				//  Average X and Y.
				double xbar = phys::maths::meanVec(x);
				double ybar = phys::maths::meanVec(y);
				//  Compute gradient and intercept
				double top = 0.0;
				double bot = 0.0;
				for ( size_t i = 0; i < n; i++ ){
					top += ( x[i] - xbar ) * ( y[i] - ybar );
					bot += ( x[i] - xbar ) * ( x[i] - xbar );
				}
				retval.first = top / bot;
				retval.second = ybar - top * xbar/ bot;
			}
			return retval;
		}

		/** Fitting data to a linear curve (read straight line) p(x) = a*x + b directly
		*/
		std::pair<double, double> linear_fit( const stdVec_d&x, const stdVec_d&y )
		{
			if(x.empty() || y.empty()) return std::pair<double,double> (0.0, 0.0);		
			assert(x.size() == y.size());
			size_t n = x.size();
			std::pair<double,double> retval;
			//  Special cases.
			if ( n == 1 ){
				retval.first = 0.0;
				retval.second = y[0];
			} else if ( n == 2 ) {
				double m = (y[1] - y[0])/(x[1] - x[0]);
				double c = y[0] - m * x[0];
				retval.first = m;
				retval.second = c;
			} else {
				stdVec_d c(4);
				for(size_t i = 0; i < n; ++i){
					c[0] += x[i];
					c[1] += x[i] * x[i];
					c[2] += y[i];
					c[3] += x[i] * y[i];
				}
				double c4 = c[0] * c[0] - c[1]*n;
			
				retval.first	= (c[0] * c[2] - c[3]*n)/c4;
				retval.second	= (c[0] * c[3] - c[1]*c[2])/c4;
			}
			return retval;
		}

		/** Fitting to a general polynomial
			We find the coefficients a0 + a1*x1 + a2*x^2 + ... + am*x^m that minimise the residuals
			aj*x^j - x_i, where x_i are our data points.
		*/
		std::pair<stdVec_d, mat> polynomial_fit( size_t m, const stdVec_d& x, const stdVec_d&y )
		{
			assert(x.size() == y.size());
			size_t n = x.size();
			mat U(m, n, 1.0);
			stdVec_d s(m), g(m), a(m), h(m);

			for(size_t i = 0; i < n; ++i){
				s[0] += U[0][i] * U[0][i];
				g[0] += x[i] * U[0][i] * U[0][i];
				a[0] += U[0][i] * y[i];
			}
			g[0] /= s[0];
			a[0] /= s[0];

			//Set up the first order polynomial U_1
			for(size_t i = 0; i < n; ++i){
				U[1][i] = x[i] * U[0][i] - g[0] * U[0][i];
				s[1] += U[1][i] * U[1][i];
				g[1] += x[i] * U[1][i] * U[1][i];
				h[1] += x[i] * U[1][i] * U[0][i];
				a[1] += U[1][i] * y[i];
			}
			g[1] /= s[1];
			h[1] /= s[0];
			a[1] /= s[1];

			//Higher order polynomials from recursive relation
			for(size_t i = 1; i < m -1; ++i){
				for(size_t j = 0; j < n; ++j){
					U[i+1][j] = x[j]*U[i][j] - g[i]*U[i][j] - h[i]*U[i-1][j];
					s[i+1] += U[i+1][j] * U[i+1][j];
					g[i+1] += x[j]*U[i+1][j]*U[i+1][j]; 
					h[i+1] += x[j]*U[i+1][j]*U[i][j];
					a[i+1] += U[i+1][j]*y[j];
				}
				g[i+1] /= s[i+1];
				h[i+1] /= s[i];
				a[i+1] /= s[i+1];
			}

			std::pair<stdVec_d,mat> retval;
			retval.first = a;
			retval.second = U;

			return retval;

		}

	}

}

