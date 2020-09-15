#include "RootSearch.h"

#include <cmath>

double func(double x) {
	//return x * x * x - 2 * x + 2;
	return cos(x) - x;
}

double deriv(double x) {
	//return 3 * x * x - 2;
	return -sin(x) - 1;
}

int main (int argc, char ** argv) {

	using NewtonRaphson = phys::roots::Newton_Raphson;

	double x0 = 0.5; //initial guess

	NewtonRaphson NR_search(&func, &deriv, x0);

	NR_search.find_root(true); //print output

	return 0;
}
