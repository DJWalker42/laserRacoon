#include "Monte_Carlo.h"

#include <cmath>
#include <iostream>



double f(double x ) {
	return sin(x); //integration: -cos(x)
}


int main (int argc, char ** argv) {

	uint32_t n {1000};

	if(argc > 1) {
		n = atoi(argv[1]);
	}

	double left {0.};
	double right {1.};

	std::pair<double, double> result {phys::MCIntegration(left, right, n,  f)};

	std::cout << "N: " << n << "\n";

	std::cout << "Analytic Soln: " << cos(0.) - cos(1.) << "\n";

	//remember, sigma provides an *estimate* of the error, which improves as n increases
	std::cout << "Soln: " << result.first << " sigma: " << result.second << std::endl;

	return 0;
}
