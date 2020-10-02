/*
 *  Program to perform an adaptive strip quadrature using the
 *  trapezoidal method.
 *
 *  If T0 defines the numerical value of the trapezium rule
 *  integration with one strip, then T1 defines it with 2
 *  strips. Generally, Tp, n = 2^p. Using the fact that the
 *  trapezium rule is of order h^2 we can obtain an estimate
 *  for the error using
 *   err = (T0 - T1)/3
 *  By comparing this estimate to the overall tolerance required
 *  we either accept the numerical integration value, or halve
 *  the strip width and repeat the process until the desired
 *  accuracy is achieved.
 *
 *  See the implementation of Adaptive::integrate and Adaptive::recursive.
 */


#include "Quadrature.h"

#include <cmath>


//normalised Lorentzian Line
/*
double f(double lambda) {

	double lambda_0 = 10.;
	double gamma = 1.;
	double x = 2 * (lambda - lambda_0) / gamma;

	return 1./(1. + x * x);
}

// analytic definite integral of the Lorentzian Line
double f_int(double x) {
	//double I0 = 1.;
	double lambda_0 = 10.;
	double gamma = 1.;

	double x_sub = 2 * (x - lambda_0) / gamma;

	return atan(x_sub) * gamma / 2.;
}
*/

double f(double x) {
	return sin(x);
}

double f_int(double x) {
	return -cos(x);
}


int main (int argc, char ** argv) {

	using Adaptive = phys::quad::Adaptive;
	using Trapezoid = phys::quad::Trapezoid;

	Adaptive adaptive;
	Trapezoid trapezoid;

	double left {0.}; // 6.
	double right {1.}; // 14.

	adaptive.set_max_count(1000);
	double adaptive_soln = adaptive.integrate(f, left, right, 1.e-3);

	std::cout << "adaptive soln: " << adaptive_soln
			<< " recursive counts: " << adaptive.get_count() << "\n";


	double soln = f_int(right) - f_int(left);

	std::cout << "number of function calls: " << adaptive.get_f_calls() << std::endl;
	std::cout << "Analytic soln: " << soln << std::endl;
	std::cout << "Adpative - soln: " << fabs(adaptive_soln - soln) << std::endl;

	for (int i = 1; i < 5; ++i) {
		uint n = static_cast<uint>(pow(10., i));
		double trapezoid_soln = trapezoid.integrate(f, left, right, n);

		std::cout << "Trap soln (" << n << ") = " << trapezoid_soln;
		std::cout << " error = " << fabs(trapezoid_soln - soln) << "\n";
	}



	return 0;
}
