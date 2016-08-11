#include <cmath>
#include <cfloat>

#include "Quadrature.h"
//Included hpp file contains the values of the knots and corresponding
//weights for both Legendre and Laguerre methods, and the defintion of the glaw array.
#include "GaussKnotsWeights.h"

namespace phys{
	namespace quad{

		/***********************************************************************************
		*	Quadrature base class
		***********************************************************************************/
		// Constructor
		Quadrature::Quadrature(uint n) : m_error(1.0), m_order(n), m_func_calls(0)
		{}

		//---------------------------------------------------------------------------------
		//	Quadrature Methods
		//---------------------------------------------------------------------------------
		double Mid_ordinate::integrate(double(*f)(double), double lft, double rht, uint N)
		{
			m_error = 0.0; 
			double h = (rht - lft)/N;
			double x = lft + h/2;
			double f1;
			double f2 = f(x); 
			double f3 = f(x+h);
			double sum = f2 + f3;
			for(uint i = 2; i < N; ++i)
			{	
				x += h;
				f1 = f2;
				f2 = f3;
				f3 = f(x + h);
				sum += f3;
				m_error += f3 - 2*f2 + f1;
			}			
			m_error *= h/24; 
			m_func_calls += N;
			return (sum*h);
		}

		double Trapeziod::integrate(double(*f)(double), double lft, double rht, uint N)
		{
			m_error = 0.0; 
			double h = (rht - lft)/N;
			double x = lft;
			double f1;
			double f2 = f(x); 
			double f3 = f(x+h);
			double sum = f2 + 2*f3; //f(lft) + 2*f(lft + h)
			for(uint i = 2; i < N; ++i)
			{	
				x += h;
				f1 = f2;
				f2 = f3;
				f3 = f(x+h);
				sum += 2*f3;
				m_error += f3 - 2*f2 + f1;			
			}
			sum += f(rht); 
			m_error *= -h/12;
			m_func_calls += N + 1;
			return (sum*h/2);
		}

		double Simpson::integrate(double(*f)(double), double lft, double rht, uint N)
		{
			double sum = 0.0;
			double h = (rht - lft)/2/N;
			double x = lft;
			double bigerr = 0.0;
			double f1, f2;
			double f3 = f(lft);
			double f4 = f(lft + h);
			double f5 = f(lft + 2*h);
			for(uint i=0; i<N; ++i)
			{
				x += 2*h;
				f1 = f3;	
				f2 = f4;	
				f3 = f5;	
				f4 = f(x+h);
				f5 = f(x+2*h);
				sum += f1 + 4*f2 + f3;
				m_error = f1 - 4*f2 + 6*f3 - 4*f4 + f5; 
				if(fabs(bigerr) < fabs(m_error)) bigerr = m_error;				
			}
			m_error = bigerr * (rht-lft)/180;
			m_func_calls += 2*N + 3;
			return (sum*h/3);
		}

		double Boole::integrate(double(*f)(double), double lft, double rht, uint N)
		{
			double sum = 0.0;
			double h = (rht - lft)/4/N;
			double x = lft;
			for(uint i = 0; i < N; ++i)
			{
				x += 4*h;
				sum += 7*f(x-4*h) + 32*f(x-3*h) + 12*f(x-2*h) 
					+ 32*f(x-h) + 7*f(x);
			}
			sum *= 2*h/45;
			m_func_calls += 5*N;
			
			if(error_wanted)
			{
				bool error_temp = error_wanted;
				error_wanted = false;
				double sum_2 = integrate(f, lft, rht, 2 * N);
				m_error = fabs(sum - sum_2)/63;
				sum = sum_2;
				error_wanted = error_temp;
			}
			return (sum);
		}


		//-------------------------------------------------------------------------------
		//	Gaussian Quadrature Methods 
		//-------------------------------------------------------------------------------
		
		/***********************************************************************************
		*	Gauss base class
		***********************************************************************************/
		//Constructor
		Gauss::Gauss(uint pts) : m_numKnots(pts), m_knots(), m_weights(), m_func_calls(0), m_precision(1.e-10)
		{}

		void Gauss::set_points(uint pts)
		{
			m_numKnots = pts;
			m_knots.clear();
			m_weights.clear();
			this->initialise();
		}

		void Gauss::set_precision(double acc)
		{
			m_precision = acc;
			m_knots.clear();
			m_weights.clear();
			this->initialise();
		}

		//Legendre ----------------------------------------------------------------------
		Legendre::Legendre(uint points) : Gauss(points) 
		{initialise();}

		void Legendre::initialise() //function called in body of constructor
		{
			// operator >> in this context is bit shift: x >> 1 means bit shift the value in x one to the right.
			uint m = (m_numKnots + 1) >> 1;
			/* Load appropriate predefined values for abscissas and weights */
			uint i = 0;
			uint max = Legen::glawsize;
			do{
				if(m_numKnots == Legen::glaw[i].n)
				{
					for(uint j = 0; j < m; ++j){
						m_knots.push_back(Legen::glaw[i].x[j]);
						m_weights.push_back(Legen::glaw[i].w[j]);
					}
				}
			}while(m_knots.empty() && ++i < max);
			/*	If the values have not been predefined for a particular n
				then compute values	*/				
			if(i == max){
				std::pair< stdVec_d, stdVec_d > xw = compute_x_w();
				m_knots = xw.first;
				m_weights = xw.second;
			}
		}

		std::pair<stdVec_d, stdVec_d> Legendre::compute_x_w()
		{
			//std::cout << "Legendre: Computing weights and abscissa for n = " << n << "\n";; 
			double x0,  x1,  dx;	/* Abscissa */
			double w0 = 0.0,  w1,  dw;	/* Weights */
			double P0, P_1, P_2;	/* Legendre polynomial values */
			double dpdx;			/* Legendre polynomial derivative */
			uint j;					/* iteration count*/
			double t0, t1, t2, t3;	/* temporary varibales to compute P values */

			/*	Add one to order and bit-shift one to right - this is so symmetry can be exploited
				in both odd and even Legendre polynomials, i.e. ignore the negative roots.
				Equivalent to:
				if(n == even) m = n/2;
				if(n == odd)  m = (n-1)/2 + 1;
			*/			
			uint m = (m_numKnots+1)>>1;
			stdVec_d vx(m), vw(m);

			/* Search for Francesco Tricomi for explanation of initial guess */
			double n2 = double(m_numKnots*m_numKnots);
			double n3 = double(n2*m_numKnots);
			t0 = 1.0 - (1.0/8.0/n2) + (1.0/8.0/n3);
			t1 = 1.0/(4.0*(double)m_numKnots+2.0);

			for (uint i = 1; i <= m; i++)
			{
				/* Find i-th root of Legendre polynomial */

				/*	Initial guess (in descending order i.e. x0[i=1] > x0[i=2] >....)
					i bit shifted 2 to the left is equivalent to multiplication by four*/
				x0 = cos(PI * double((i<<2)-1) * t1) * t0;

				/* Newton iterations, at least one */
				j = 0;
				dx = dw = DBL_MAX;
				do 
				{
					/* Compute Legendre polynomial value at x0 */
					P_1 = 1.0;
					P0  = x0;
#if 0
					/* Simple, not optimized version */
					for (uint k = 2; k <= m_numKnots; k++)
					{
						P_2 = P_1;
						P_1 = P0;
						t2 = x0*P_1;
						t3 = (double)(k-1)/(double)k;
						P0 = t2 + t3*(t2 - P_2);
					}
#else
					/* Optimized version using lookup tables */
					if (m_numKnots<1024)
					{
						/* Use fast algorithm for small n*/
						for (uint k = 2; k <= m_numKnots; k++)
						{
							P_2 = P_1;
							P_1 = P0;
							t2  = x0*P_1;
							P0 = t2 + Legen::ltbl[k]*(t2 - P_2);
						}
					}else{
						/* Use general algorithm for other n */
						for (uint k = 2; k < 1024; k++)
						{
							P_2 = P_1;
							P_1 = P0;
							t2  = x0*P_1;
							P0 = t2 + Legen::ltbl[k]*(t2 - P_2);
						}
						for (uint k = 1024; k <= m_numKnots; k++)
						{
							P_2 = P_1;
							P_1 = P0;
							t2 = x0*P_1;
							t3 = (double)(k-1)/(double)k;					
							P0 = t2 + t3*(t2 - P_2);
						}
					}
#endif
					/* Compute Legendre polynomial derivative at x0 */
					dpdx = ((x0 * P0 - P_1) * double(m_numKnots))/(x0*x0-1.0);

					/* Newton step */
					x1 = x0 - P0 / dpdx;

					/* Weight computing */
					w1 = 2.0 / ((1.0 - x1 * x1) * dpdx * dpdx);

					/* Compute weight w0 on first iteration, needed for dw */
					if (j == 0) w0 = 2.0 / ((1.0 - x0 * x0) * dpdx * dpdx);

					dx = x0 - x1;
					dw = w0 - w1;

					x0 = x1;
					w0 = w1;
					j++;
				} while((fabs(dx) > m_precision || fabs(dw) > m_precision) && j < 100);

				/*	Remember that the roots are found in descending order but the 
					integration is expecting values in ascending order.	*/
				vx[(m-1) - (i-1)] = x1;
				vw[(m-1) - (i-1)] = w1;
			}//end of outer for loop
			std::pair<stdVec_d, stdVec_d> retval(vx, vw);
			return retval;
		}

		double Legendre::integrate(double(*f)(double), double lft, double rht)
		{			
			double A = 0.5 * (rht - lft);
			double B = 0.5 * (rht + lft);
			/*	Bit shift (n+1) one to the right 
				-- exploits symmetry of Legendre polynomials for both odd and even n */ 
			uint m = (m_numKnots + 1)>>1;
			m_func_calls += 2 * m;
			double s, Ax;
			if(m_numKnots & 1) /* n - odd */
			{
				--m_func_calls;
				s = m_weights[0]*(f(B));
				for (uint i = 1; i < m; i++)
				{
					Ax = A * m_knots[i];
					s += m_weights[i] * (f(B + Ax) + f(B - Ax));
				}
			}else{ /* n - even */	
				s = 0.0;
				for (uint i = 0; i < m; i++)
				{
					Ax = A * m_knots[i];
					s += m_weights[i] * (f(B + Ax) + f(B - Ax));			
				}
			}
			return A * s;
		}

		//Laguerre -----------------------------------------------------------------
		Laguerre::Laguerre( uint points, bool modify ) : Gauss(points), modified(modify) 
		{initialise();}

		void Laguerre::initialise()//function called in body of constructor
		{
			/* Load appropriate predefined values for abscissa and weights */
			uint i = 0;
			uint max = Lague::glawsize;
			do{
				if (m_numKnots == Lague::glaw[i].n)
				{
					for(uint j = 0; j < m_numKnots; ++j){
						m_knots.push_back(Lague::glaw[i].x[j]);
						m_weights.push_back(Lague::glaw[i].w[j]);
					}
				}
			}while(m_knots.empty() && ++i < max);
			/*	If the values have not been predefined for a particular n
				then compute values	*/				
			if(i == max){	
				std::pair<stdVec_d, stdVec_d> xw = compute_x_w();
				m_knots = xw.first;
				m_weights = xw.second;
			}

			if(modified)
				for(uint i = 0; i < m_knots.size(); ++i)
					m_weights[i] *= exp(m_knots[i]);
		}

		std::pair<stdVec_d, stdVec_d> Laguerre::compute_x_w()
		{	
			//std::cout << "Laguerre: Computing weights and abscissa for n = " << m_numKnots << "\n"; 
			double x0, x1, dx;
			double w0 = 0.0, w1, dw;
			double L_0, L_1, L_2;
			double dldx;
			double t1, t2, t3;
			uint j;

			stdVec_d vx(m_numKnots), vw(m_numKnots);

			for (uint i = 1; i <= m_numKnots; i++) {
				if (i == 1)       // initial approximation for zeros (Stroud & Secrest)
					x0 = 3./(1. + 2.4 * double(m_numKnots));             // 1st zero
				else if (i == 2)
					x0 = 15./(1. + 2.5 * double(m_numKnots)) + vx[0];    // 2nd zero
				else {
					L_0 = (1./double(i) + 2.55)/1.9;			// recurrence
					x0 = (1. + L_0) * vx[i-2] - L_0*vx[i-3];
				}
				j = 0;
				dx = dw = DBL_MAX;

				do{
					L_0 = 1.0 - x0;
					L_1 = 1.0;
					for(uint k = 2; k<=m_numKnots; ++k){
						L_2 = L_1;
						L_1 = L_0;
						t1 = 2*L_1 - L_2;
						t2 = L_1*(1.0 + x0);
						t3 = 1.0/double(k); 
						L_0 = t1 + t3*(L_2 - t2);
					}
					dldx = x0 ? double(m_numKnots)*(L_0-L_1)/x0 : -double(m_numKnots)*L_0; //avoids divide by zero

					x1 = x0 - L_0/dldx;
					w1 = 1/x1/dldx/dldx;
					if(j==0) w0 = 1/x0/dldx/dldx;

					dx = x0-x1;
					dw = w0-w1;
					
					x0 = x1;
					w0 = w1;

				}while( (fabs(dx)>m_precision || fabs(dw) > m_precision) && ++j <100);

				vx[i-1] = x1;
				vw[i-1] = w1;
			}

			std::pair<stdVec_d, stdVec_d> retval(vx, vw);
			return retval;
		}

		double Laguerre::integrate( double(*f)(double), double lft, double )
		{			
			double sum = 0.0;	
			for(uint i = 0; i < m_numKnots; ++i)
				sum += m_weights[i] * f(lft+m_knots[i]);	

			m_func_calls += m_numKnots;
			return sum;
		}

		/*******************************************************************************************************
		*	Other methods
		*******************************************************************************************************/
		//Adaptive extension --------------------------------------------------------------------------------
		double Adaptive::integrate(double(*f)(double), double lft, double rht, uint pts, double tol)
		{
			if (tol != double()) m_tolerance = tol;
			m_I = f;
			reset_count();
			double m_error = 0.0;
			double retval;

			retval = recursive(lft, rht, pts, m_tolerance, m_count, m_error);
			m_num_f_calls = m_pQuad->get_func_calls();

			if (m_count > m_cnt_max)
				std::cout << "Warning: max count reached.\n";
			if (m_print_output)
				std::cout << "Solution: " << retval << " +/- " << m_error << std::endl;

			return retval;
		}

		double Adaptive::recursive(double lft, double rht, uint strips, double tol, uint&cnt, double&m_error)
		{
			double T0 = m_pQuad->integrate(m_I, lft, rht, strips);
			double T1 = m_pQuad->integrate(m_I, lft, rht, strips * 2);

			double mid = (rht + lft) / 2;
			double sum = 0.0;

			double err = fabs(T1 - T0) / (pow(2, m_order) - 1.0);
			if ((err <= tol && cnt > m_cnt_min) || cnt > m_cnt_max) // accept the integration
			{
				sum = T1;
				m_error += err;
				return sum;
			}
			else
			{
				sum += recursive(lft, mid, strips, tol/2., ++cnt, m_error);
				sum += recursive(mid, rht, strips, tol/2., cnt, m_error);
				return sum;
			}

		}

		//Romberg integration ----------------------------------------------------------------------------------
		double Romberg::integrate(double(*f)(double), double lft, double rht, double tolerance, uint max_level)
		{
			if (tolerance != double()) m_tol = tolerance;
			if (max_level != uint()) m_level = max_level;

			unsigned lvl = m_level, k = 0;
			stdVec_d s(lvl, 1.0);
			double temp = 0.0;

			Trapeziod trap;

			m_error = 1.0;
			while (++k < s.size() && m_error > m_tol)
			{
				for (unsigned i = 1; i <= k; i++)
				{
					if (i == 1)
					{
						temp = s[i];
						s[i] = trap.integrate(f, lft, rht, uint(pow(2, k - 1)));
					}
					else
					{
						s[k] = (pow(4, i - 1)*s[i - 1] - temp) / (pow(4, i - 1) - 1);
						temp = s[i];
						s[i] = s[k];
					}
				}
				if (k > 1){
					m_error = fabs(s[k] - s[k - 1]);
					//force another extrapolation level if m_error is suspiciously small
					if (m_error < 1000 * DBL_EPSILON) m_error = 1.0;
				}
			}
			if (k == lvl){
				std::cout << "Warning: Max extrapolation level reached.\n";
				m_error = fabs(s[k - 1] - s[k - 2]);
			}
			if (m_print_output)
				std::cout << "Solution: " << s[k - 1] << " +/- " << m_error
				<< "in " << k << " levels." << std::endl;

			m_func_calls += trap.get_func_calls();

			return s[k - 1];
		}
	}
}
