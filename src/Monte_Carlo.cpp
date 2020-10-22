#include "Monte_Carlo.h"

#include "RandomNumGen.h"
#include "Maths.h"

namespace phys {

std::pair<double, double> MCIntegration(double a, double b, uint32_t n, double(*f)(double)) {

	assert(n > 0);
	assert(a < b);

	RNG<> generator(a,b);

	double sum {0.};
	double sum2 {0.};

	for (int i = 0; i < n; ++i) {
		double val = f(generator.random_number());
		sum += val;
		sum2 += val * val;
	}

	double s {sum/double(n)};
	double s2 {sum2/double(n)};

	double sigma {sqrt( (s2 - (s*s)) / double(n - 1) )};

	return { s * (b - a), sigma * (b - a) };
}



}//namespace
