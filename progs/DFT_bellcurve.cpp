#include "Fourier.h"

#include <cmath>

double bellcurve(double x) {
	double mu {.5};
	double si {.1};
	double pi {4. * atan(1.)};

	return exp(-(x-mu)*(x-mu)/2./si/si)/si/sqrt(2.*pi);
}

int main (int argc, char ** argv) {

	const int n {64};

	double f0 {1./n};

	double step {1./(n - 1)};

	phys::stdVec_c f(n); //vector of complex (elements initialised to zero)

	double x {0.};

	for(auto& elem : f) {
		elem = bellcurve(x); //assigns to real part of complex
		x += step;
	}

	//Fourier transform of f
	phys::stdVec_c g(phys::DFT(f)); //uses copy constructor

	//normalise the transform, and take complex conjugate
	for (auto& elem : g) {
		elem *= f0;
		elem = conj(elem);
	}

	//perform the inverse transform
	phys::stdVec_c ft(phys::DFT(g));

	//ft should equal f
	fprintf(stdout, "#f\tft\n");
	for (int i = 0; i < n; ++i) {
		fprintf(stdout, "%f\t%f\n", f[i].real(), ft[i].real());
	}

	return 0;
}
