#include "Quadrature.h"
#include "FormatOutput.h"

int main (int argc, char ** argv) {

	double analytic = exp(4.) - exp(1.);

	printf("n\terror\n");

	for (int n = 2; n < 9; ++n) {
		phys::Legendre leg(n);
		double result = leg.integrate(exp, 1., 4.);

		printf("%d\t%s\n", n, phys::format(fabs(analytic - result), 3).c_str());
	}

	return 0;
}
