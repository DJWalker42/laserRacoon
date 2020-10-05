#include "GaussSeidel.h"

namespace phys {

void gaussSeidel(	stdVec_d& f,
					const std::array<double,4> coeffs,
					double x_1,
					double x_n,
					double tol,
					int max_iter ) {

	double a {coeffs[0]}, b {coeffs[1]}, c {coeffs[2]}, d {coeffs[3]};

	int n {int(f.size())};
	double h {(x_n - x_1)/(n - 1)};

	double theta {c - 2 * a /h/h};
	double phi {(a/h/h) + (b/2/h)};
	double psi {(a/h/h) - (b/2/h)};

	int iter_count{0};

	bool tol_achieved {false};

	double x {0.};
	double ff {0.};

	while ( tol_achieved == false && ++iter_count < max_iter) {

		tol_achieved = true;

		for (int i = 1; i < n - 1; ++i) {
			x = x_1 + h * i;
			ff = -(phi * f[i+1] + psi * f[i-1] - d * x) / theta;

			if ( fabs( (ff - f[i])/ff ) > tol ) tol_achieved = false;

			f[i] = ff;
		}
	}

	/*
	 * TODO: add warning message here should iter_count reach max_iter
	 */
}

} //namespace
