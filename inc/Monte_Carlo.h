#ifndef MONTE_CARLO
#define MONTE_CARLO

#include <utility>
#include <cstdint>

namespace phys {


//1D monte carlo integration of function f over [a,b] using n randomly sampled points
//returns the integration result and estimate of the error
std::pair<double, double> MCIntegration(double a, double b, uint32_t n, double(*f)(double));


} //namespace

#endif //header guard
