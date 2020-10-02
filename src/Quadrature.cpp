#include <cmath>
#include <cfloat>

#include "Quadrature.h"
//Included hpp file contains the values of the knots and corresponding
//weights for both Legendre and Laguerre methods, and the definition of the glaw array.
#include "GaussKnotsWeights.h"

namespace phys {
namespace quad {


/*
 * |E| is the upper bound of the error of the quadrature rule
 * M_p is the upper bound of pth derivative in [a,b]
 * N is the number of "strips"
 *
 * Midpoint:
 * 	|E| <= M_2 * (b - a)^3 / 24 / N^2
 *
 * Trapezoid:
 *  |E| <= M_2 * (b - a)^3 / 12 / N^2
 *
 * Simpson:
 * 	|E| <= M_4 * (b - a)^5 / 180 / N^4
 *
 * 	Boole:
 * 	 ???
 *
 * Finite difference (central):
 *  f'2	~ (f(x+h) - 2*f(x) + f(x-h)) / h^2
 *  f'4 ~ (f(x+2h) - 4*f(x+h) + 6*f(x) - 4*f(x-h) + f(x-2h)) / h^4
 *  f'6 ~ ???
 */


/***********************************************************************************
 *	Quadrature base class
 ***********************************************************************************/
// Constructor
Quadrature::Quadrature(uint n) :
		m_error(1.0), m_order(n), m_func_calls(0) {
}

//---------------------------------------------------------------------------------
//	Quadrature Methods
//---------------------------------------------------------------------------------
double MidOrdinate::integrate(double (*f)(double), double a, double b) {
	double h { b - a };
	double fmid { f((b + a) / 2.) };

	//taken to estimate the error (may want to put this elsewhere)
	double fa { f(a) };
	double fb { f(b) };
	//finite diff. 2nd deriv, are the step lengths the same?
	m_error = -(fb - 2 * fmid + fa) * h / 24.;

	m_func_calls += 1;

	return fmid * h;
}

double MidOrdinate::integrate(double (*f)(double), double a, double b, uint N) {

	if (N == 0) throw std::runtime_error("MidOrdinate::integrate: N cannot be zero");

	if (N == 1) {
		return integrate(f, a, b);
	} else {
		//N > 1
		m_error = 0.0;
		double h {(b - a) / N};
		double x {a + h / 2};
		double f1;
		double f2 {f(x)};
		x += h;
		double f3 {f(x)};
		double sum {f2 + f3};
		double bigerr {0.0};
		for (uint i = 2; i < N; ++i) {
			x += h;
			f1 = f2;
			f2 = f3;
			f3 = f(x);
			sum += f3;
			//finite diff. 2nd deriv
			m_error += f3 - 2 * f2 + f1;
			if (fabs(bigerr) < fabs(m_error))
				bigerr = m_error;	 //estimate of the error upper-bound on the interval [a,b]
		}
		m_error = bigerr * -(b - a) / 24; //check this
		m_func_calls += N;
		return sum * h;
	}
}

double Trapezoid::integrate(double(*f)(double), double a, double b) {
	double h {b - a};

	//function evaluations at extremes of interval
	double fa {f(a)};
	double fb {f(b)};

	//midpoint for error estimate
	double fmid {f((a + b) / 2.)};

	//finite difference 2nd derivative (error order 2)
	m_error = -(fb - 2 * fmid + fa) * h / 12.; //is this correct? check step length of the FD

	m_func_calls += 2;
	//single strip estimate of the integral
	return (fa + fb) * h / 2;
}

double Trapezoid::integrate(double (*f)(double), double a, double b, uint N) {

	if (N == 0) throw std::runtime_error("Trapezoid::integrate: N cannot be zero");

	if (N == 1) {
		return integrate(f, a, b);
	} else {
		//N > 1
		m_error = 0.0;
		double h {(b - a) / N};
		double x {a};
		double f1;
		double f2 {f(x)};
		x += h;
		double f3 {f(x)};
		double sum {f2 + 2 * f3};
		double bigerr {0.0};
		for (uint i = 2; i < N; ++i) {
			x += h;
			f1 = f2;
			f2 = f3;
			f3 = f(x);
			sum += 2 * f3;
			//using finite difference 2nd derivative
			m_error += f3 - 2 * f2 + f1;
			if (fabs(bigerr) < fabs(m_error)) bigerr = m_error;
		}
		sum += f(b);
		m_error = bigerr * (b - a) / N/ 12; //check this is correct
		m_func_calls += N + 1;
		return sum * h / 2;
	}
}

double Simpson::integrate(double(*f)(double), double a, double b) {
	//Here we could just call the composite rule with N = 1
	double h {b - a};
	double c {(a + b) / 2.};

	//to estimate the error we need a finite difference scheme for the 4th derivative

	m_func_calls += 3;

	return h * (f(a) + 4 * f(c) + f(b)) / 6.;

}

double Simpson::integrate(double (*f)(double), double a, double b, uint N) {
	double sum {0.0};
	double h {(b - a) / 2 / N}; //we've halved the usual definition of h
	double x {a};
	double bigerr {0.0};
	double f1, f2;
	double f3 {f(a)};
	double f4 {f(a + h)};
	double f5 {f(a + 2 * h)};
	for (uint i = 0; i < N; ++i) {
		x += 2 * h;
		f1 = f3;
		f2 = f4;
		f3 = f5;
		f4 = f(x + h);
		f5 = f(x + 2 * h);
		sum += f1 + 4 * f2 + f3;
		//finite diff. 4th deriv 2nd order accurate
		m_error = f1 - 4 * f2 + 6 * f3 - 4 * f4 + f5;
		if (fabs(bigerr) < fabs(m_error)) bigerr = m_error;
	}
	m_error = bigerr * -(b - a) / 180;
	m_func_calls = 2 * N + 3;
	return (sum * h / 3);
}

double Boole::integrate(double(*f)(double), double a, double b) {
	//instead of computing the primitive integral we just call composite with N = 1
	return integrate(f, a, b, 1);
}

double Boole::integrate(double (*f)(double), double a, double b, uint N) {
	double sum {0.0};
	double h {(b - a) / 4 / N};
	double x {a};
	for (uint i = 0; i < N; ++i) {
		x += 4 * h;
		sum += 7 * f(x - 4 * h) + 32 * f(x - 3 * h) + 12 * f(x - 2 * h)
				+ 32 * f(x - h) + 7 * f(x);
	}
	sum *= 2 * h / 45;
	m_func_calls += 5 * N;

	//here we estimate the error by performing the integration with twice the
	//number of strips (h/2). Instead, what order of derivative do we need to
	//approximate to estimate the error?

	if (error_wanted) {
		bool error_temp {error_wanted};
		error_wanted = false;
		double sum_2 = integrate(f, a, b, 2 * N);
		m_error = fabs(sum - sum_2) / 63; //Boole's rule is order 6 accurate
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
Gauss::Gauss(uint pts) :
		m_numKnots(pts), m_knots(), m_weights(),
		m_func_calls(0), m_precision(1.e-10) {}

void Gauss::set_points(uint pts) {
	m_numKnots = pts;
	m_knots.clear();
	m_weights.clear();
	this->initialise();
}

void Gauss::set_precision(double acc) {
	m_precision = acc;
	m_knots.clear();
	m_weights.clear();
	this->initialise();
}

//Legendre ----------------------------------------------------------------------
Legendre::Legendre(uint points) :
		Gauss(points) {
	initialise();
}

void Legendre::initialise() //function called in body of constructor
{
	// operator >> in this context is bit shift: x >> 1 means bit shift the value in x one to the right.
	uint m = (m_numKnots + 1) >> 1;
	/* Load appropriate predefined values for abscissas and weights */
	uint i = 0;
	uint max = Legen::glawsize;
	do {
		if (m_numKnots == Legen::glaw[i].n) {
			for (uint j = 0; j < m; ++j) {
				m_knots.push_back(Legen::glaw[i].x[j]);
				m_weights.push_back(Legen::glaw[i].w[j]);
			}
		}
	} while (m_knots.empty() && ++i < max);
	/*	If the values have not been predefined for a particular n
	 then compute values	*/
	if (i == max) {
		std::pair<stdVec_d, stdVec_d> xw = compute_x_w();
		m_knots = xw.first;
		m_weights = xw.second;
	}
}

std::pair<stdVec_d, stdVec_d> Legendre::compute_x_w() {
	//std::cout << "Legendre: Computing weights and abscissa for n = " << n << "\n";;
	double x0, x1, dx; /* Abscissa */
	double w0 = 0.0, w1, dw; /* Weights */
	double P0, P_1, P_2; /* Legendre polynomial values */
	double dpdx; /* Legendre polynomial derivative */
	uint j; /* iteration count*/
	double t0, t1, t2, t3; /* temporary variables to compute P values */

	/*	Add one to order and bit-shift one to right - this is so symmetry can be exploited
	 in both odd and even Legendre polynomials, i.e. ignore the negative roots.
	 Equivalent to:
	 if(n == even) m = n/2;
	 if(n == odd)  m = (n-1)/2 + 1;
	 */
	uint m = (m_numKnots + 1) >> 1;
	stdVec_d vx(m), vw(m);

	/* Search for Francesco Tricomi for explanation of initial guess */
	double n2 = double(m_numKnots * m_numKnots);
	double n3 = double(n2 * m_numKnots);
	t0 = 1.0 - (1.0 / 8.0 / n2) + (1.0 / 8.0 / n3);
	t1 = 1.0 / (4.0 * (double) m_numKnots + 2.0);

	for (uint i = 1; i <= m; i++) {
		/* Find i-th root of Legendre polynomial */

		/*	Initial guess (in descending order i.e. x0[i=1] > x0[i=2] >....)
		 i bit shifted 2 to the left is equivalent to multiplication by four*/
		x0 = cos(PI * double((i << 2) - 1) * t1) * t0;

		/* Newton iterations, at least one */
		j = 0;
		dx = dw = DBL_MAX;
		do {
			/* Compute Legendre polynomial value at x0 */
			P_1 = 1.0;
			P0 = x0;
#if 0
					/* Simple, non optimised version */
					for (uint k = 2; k <= m_numKnots; k++)
					{
						P_2 = P_1;
						P_1 = P0;
						t2 = x0*P_1;
						t3 = (double)(k-1)/(double)k;
						P0 = t2 + t3*(t2 - P_2);
					}
#else
			/* Optimised version using lookup tables */
			if (m_numKnots < 1024) {
				/* Use fast algorithm for small n*/
				for (uint k = 2; k <= m_numKnots; k++) {
					P_2 = P_1;
					P_1 = P0;
					t2 = x0 * P_1;
					P0 = t2 + Legen::ltbl[k] * (t2 - P_2);
				}
			} else {
				/* Use general algorithm for other n */
				for (uint k = 2; k < 1024; k++) {
					P_2 = P_1;
					P_1 = P0;
					t2 = x0 * P_1;
					P0 = t2 + Legen::ltbl[k] * (t2 - P_2);
				}
				for (uint k = 1024; k <= m_numKnots; k++) {
					P_2 = P_1;
					P_1 = P0;
					t2 = x0 * P_1;
					t3 = (double) (k - 1) / (double) k;
					P0 = t2 + t3 * (t2 - P_2);
				}
			}
#endif
			/* Compute Legendre polynomial derivative at x0 */
			dpdx = ((x0 * P0 - P_1) * double(m_numKnots)) / (x0 * x0 - 1.0);

			/* Newton step */
			x1 = x0 - P0 / dpdx;

			/* Weight computing */
			w1 = 2.0 / ((1.0 - x1 * x1) * dpdx * dpdx);

			/* Compute weight w0 on first iteration, needed for dw */
			if (j == 0)
				w0 = 2.0 / ((1.0 - x0 * x0) * dpdx * dpdx);

			dx = x0 - x1;
			dw = w0 - w1;

			x0 = x1;
			w0 = w1;
			j++;
		} while ((fabs(dx) > m_precision || fabs(dw) > m_precision) && j < 100);

		/*	Remember that the roots are found in descending order but the
		 integration is expecting values in ascending order.	*/
		vx[(m - 1) - (i - 1)] = x1;
		vw[(m - 1) - (i - 1)] = w1;
	} //end of outer for loop
	std::pair<stdVec_d, stdVec_d> retval(vx, vw);
	return retval;
}

double Legendre::integrate(double (*f)(double), double lft, double rht) {
	double A = 0.5 * (rht - lft);
	double B = 0.5 * (rht + lft);
	/*	Bit shift (n+1) one to the right
	 -- exploits symmetry of Legendre polynomials for both odd and even n */
	uint m = (m_numKnots + 1) >> 1;
	m_func_calls += 2 * m;
	double s, Ax;
	if (m_numKnots & 1) /* n - odd */
	{
		--m_func_calls;
		s = m_weights[0] * (f(B));
		for (uint i = 1; i < m; i++) {
			Ax = A * m_knots[i];
			s += m_weights[i] * (f(B + Ax) + f(B - Ax));
		}
	} else { /* n - even */
		s = 0.0;
		for (uint i = 0; i < m; i++) {
			Ax = A * m_knots[i];
			s += m_weights[i] * (f(B + Ax) + f(B - Ax));
		}
	}
	return A * s;
}

//Laguerre -----------------------------------------------------------------
Laguerre::Laguerre(uint points, bool modify) :
		Gauss(points), modified(modify) {
	initialise();
}

void Laguerre::initialise() //function called in body of constructor
{
	/* Load appropriate predefined values for abscissa and weights */
	uint i = 0;
	uint max = Lague::glawsize;
	do {
		if (m_numKnots == Lague::glaw[i].n) {
			for (uint j = 0; j < m_numKnots; ++j) {
				m_knots.push_back(Lague::glaw[i].x[j]);
				m_weights.push_back(Lague::glaw[i].w[j]);
			}
		}
	} while (m_knots.empty() && ++i < max);
	/*	If the values have not been predefined for a particular n
	 then compute values	*/
	if (i == max) {
		std::pair<stdVec_d, stdVec_d> xw = compute_x_w();
		m_knots = xw.first;
		m_weights = xw.second;
	}

	if (modified)
		for (uint i = 0; i < m_knots.size(); ++i)
			m_weights[i] *= exp(m_knots[i]);
}

std::pair<stdVec_d, stdVec_d> Laguerre::compute_x_w() {
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
			x0 = 3. / (1. + 2.4 * double(m_numKnots));             // 1st zero
		else if (i == 2)
			x0 = 15. / (1. + 2.5 * double(m_numKnots)) + vx[0];    // 2nd zero
		else {
			L_0 = (1. / double(i) + 2.55) / 1.9;			// recurrence
			x0 = (1. + L_0) * vx[i - 2] - L_0 * vx[i - 3];
		}
		j = 0;
		dx = dw = DBL_MAX;

		do {
			L_0 = 1.0 - x0;
			L_1 = 1.0;
			for (uint k = 2; k <= m_numKnots; ++k) {
				L_2 = L_1;
				L_1 = L_0;
				t1 = 2 * L_1 - L_2;
				t2 = L_1 * (1.0 + x0);
				t3 = 1.0 / double(k);
				L_0 = t1 + t3 * (L_2 - t2);
			}
			dldx = x0 ?
					double(m_numKnots) * (L_0 - L_1) / x0 :
					-double(m_numKnots) * L_0; //avoids divide by zero

			x1 = x0 - L_0 / dldx;
			w1 = 1 / x1 / dldx / dldx;
			if (j == 0)
				w0 = 1 / x0 / dldx / dldx;

			dx = x0 - x1;
			dw = w0 - w1;

			x0 = x1;
			w0 = w1;

		} while ((fabs(dx) > m_precision || fabs(dw) > m_precision) && ++j < 100);

		vx[i - 1] = x1;
		vw[i - 1] = w1;
	}

	std::pair<stdVec_d, stdVec_d> retval(vx, vw);
	return retval;
}

double Laguerre::integrate(double (*f)(double), double lft, double) {
	double sum = 0.0;
	for (uint i = 0; i < m_numKnots; ++i)
		sum += m_weights[i] * f(lft + m_knots[i]);

	m_func_calls += m_numKnots;
	return sum;
}

/*******************************************************************************************************
 *	Other methods
 *******************************************************************************************************/
//Adaptive extension --------------------------------------------------------------------------------
double Adaptive::integrate(double (*f)(double), double a, double b, double tol) {

	m_I = f;
	reset_count();
	double retval;

	double T0 {m_pQuad->integrate(f, a, b)};

	retval = recursive(a, b, tol, T0);

	if (m_count > m_cnt_max)
		std::cout << "Warning: max count reached.\n";

	return retval;
}

double Adaptive::recursive(double a, double b, double tol, double whole) {

	m_count++;

	double mid = (a + b) / 2.;

	double left {m_pQuad->integrate(m_I, a, mid)};
	double right {m_pQuad->integrate(m_I, mid, b)};

	double delta {left + right - whole};
	double order {pow(2., m_order) - 1.0};

	//condition to terminate recursion
	if (fabs(delta) <= order * tol || m_count > m_cnt_max) {
		return left + right + delta / order; //add in the error correction
	}

	return recursive(a, mid, tol/2., left) +
			recursive(mid, b, tol/2., right);
}

//Romberg integration ----------------------------------------------------------------------------------
double Romberg::integrate(double (*f)(double), double lft, double rht,
		double tolerance, uint max_level) {
	if (tolerance != double())
		m_tol = tolerance;
	if (max_level != uint())
		m_level = max_level;

	unsigned lvl = m_level, k = 0;
	stdVec_d s(lvl, 1.0);
	double temp = 0.0;

	Trapezoid trap;

	m_error = 1.0;
	while (++k < s.size() && m_error > m_tol) {
		for (unsigned i = 1; i <= k; i++) {
			if (i == 1) {
				temp = s[i];
				s[i] = trap.integrate(f, lft, rht, uint(pow(2, k - 1)));
			} else {
				s[k] = (pow(4, i - 1) * s[i - 1] - temp) / (pow(4, i - 1) - 1);
				temp = s[i];
				s[i] = s[k];
			}
		}
		if (k > 1) {
			m_error = fabs(s[k] - s[k - 1]);
			//force another extrapolation level if m_error is suspiciously small
			if (m_error < 1000 * DBL_EPSILON)
				m_error = 1.0;
		}
	}
	if (k == lvl) {
		std::cout << "Warning: Max extrapolation level reached.\n";
		m_error = fabs(s[k - 1] - s[k - 2]);
	}
	if (m_print_output)
		std::cout << "Solution: " << s[k - 1] << " +/- " << m_error << "in "
				<< k << " levels." << std::endl;

	m_func_calls += trap.get_func_calls();

	return s[k - 1];
}

}//namespace quad
}//namespace phys
