#ifndef LASER_RACOON_GAUSS_SEIDEL_H
#define LASER_RACOON_GAUSS_SEIDEL_H

#include "DynVector.h"
#include <array>

namespace phys {

void gaussSeidel(	stdVec_d& f,  //input: initial guess of solution, output: improved solution
					const std::array<double,4> coeffs,  //general 2nd order ODE coefficients
					double x_start, // start of domain
					double x_end, // end of domain
					double tol, // error tolerance in iterated solution
					int max_iter //maximum iterations to perform
);



} //namespace

#endif //header guard
