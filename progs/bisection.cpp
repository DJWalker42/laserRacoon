#include "RootSearch.h"

#include <cmath>

double func(double x) {
	return cos(x) - x;
}

int main (int argc, char ** argv) {

	//we know a single root to the function cos(x) - x lies between 0 and pi/2

	// there is a bug in RootSearch::check_brackets_set() that prevents us from
	// using zero as a valid initial bracket value, you should fix it.
	double left {0.1};
	double right {1.57};

	phys::roots::Bisection bisection(&func, left, right);

	bisection.find_root(true); //we want to print the output

	return 0;
}
