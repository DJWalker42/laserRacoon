#include "RandomNumGen.h"
#include "Visualise.h"

/*
 *  Monte-Carlo dart board experiment.
 *
 *  1000 darts are randomly "thrown" at a unit square board. An estimate for
 *  pi can be made by counting the number of darts that land within the quarter
 *  unit circle, centred on bottom left corner of the square.
 *
 *  Here we store data every 10 throws such that it can be plotted to show convergence
 *  (or lack thereof) on to the actual value of pi.
 */

int main (int argc, char ** argv) {

	const int n {1000}; //number of darts to throw
	const int c {10}; //print data every c darts thrown

	phys::RNG<> rng(0., 1.); //uniform distribution on [0, 1)

	int hits {0};

	phys::stdVec_d darts_thrown;
	phys::stdVec_d pi_estimate;

	for (int i = 0; i < n; ++i) {
		double x {rng.random_number()};
		double y {rng.random_number()};

		double r = x * x + y * y; //unit circle so no sqrt required

		if ( r <= 1 ) ++hits;

		if (i != 0 && i % c == 0) {
			darts_thrown.push_back(i);
			pi_estimate.push_back(4. * hits / i);
		}
	}

	printf("\nfinal pi estimate: %f\n\n", 4.* hits/n);

	phys::visual::Viewer viewer;
	viewer.set_y_range(2.5, 3.5);
	//viewer.withLines();
	viewer.plot(darts_thrown, pi_estimate);
	//hit any key to see the pi line
	viewer.draw_line(0., 4*atan(1.));
}
