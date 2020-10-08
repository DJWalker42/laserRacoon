#ifndef DATAFIT_HPP
#define DATAFIT_HPP

#include "DynVector.h"
#include "DynMatrix.h"

namespace phys {

/**	Linear Least Squares fit to data supplied by the user; assumes data is linear.
 @param x vector<double> independent data points
 @param y vector<double> corresponding dependent data points
 @returns a pair of doubles; first is the gradient computed, second is the y-axis intercept
 */
std::pair<double, double> llsq(const stdVec_d &x, const stdVec_d &y);

/**	Direct linear fit to the data supplied by the user
 @param x vector<double> independent data points
 @param y vector<double> corresponding dependent data points
 @returns a pair of doubles; first is the gradient computed, second is the y-axis intercept
 */
std::pair<double, double> linear_fit(const stdVec_d &x, const stdVec_d &y);

/**	Polynomial fit to the data supplied by the user
 @param m (max) polynomial order to which to fit the data (e.g. m = 2 means fit data with a quadratric)
 @param x vector<double> independent data points
 @param y vector<double> corresponding dependent data points
 @returns a pair;	first is a vector<double> of the polynomial coefficients (a_m).
 second is a matrix<double>, dimensions (m,n) of the generated polynomials
 */
std::pair<stdVec_d, mat_d> polynomial_fit(size_t m, const stdVec_d &x, const stdVec_d &y);

} //namespace

#endif //header guard
