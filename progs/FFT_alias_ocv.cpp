#include "Fourier.h"
#include "Visualise.h"
#include "error_functions.h"

#include <cmath>

//doesn't exist on Windows machines
#include <unistd.h>  //getopt, optarg, optopt

struct options_t {
	double time_range;
};

//an option letter followed by a colon means option requires a value
#define OPTION_FMT "r:h"
#define USAGE_FMT "%s [-r <time-range>]"


double f(double x) {
	return cos(5. * M_PI * x);
}

/*
 * Number of data points (samples) is 32. To increase your sample rate
 * you decrease the time range.
 *
 * Note that:
 * time ranges that don't precisely give a discrete frequency of 2.5Hz cause "leakage"
 * given the number of data points used at what time range does aliasing occur?
 *
 */

int main(int argc, char ** argv) {

	options_t options;
	options.time_range = 1.;

	int c;
	while ((c = getopt(argc, argv, OPTION_FMT)) != -1) {
	    switch (c) {
	    case 'r':
	    	options.time_range = atof(optarg);
	    	break;
	    case '?':
	        if (optopt == 'r') {
	            phys::errExit("Option -f requires a value");
	        } else if ( isprint(optopt) ) {
	            phys::errExit("Option '-%c' unknown", optopt);
	        } else {
	            phys::errExit("Option character '\\x%x' unknown", optopt);
	        }
	        break;
	    case 'h': /* fall through */
	    default:
	    	phys::usageErr(USAGE_FMT, argv[0]); //terminates
	    	break;
	    }
	}

	const int n {32};

	phys::stdVec_c data(n);
	phys::stdVec_d freq_axis;

	double time_step = options.time_range / double(n);

	for (int i = 0; i < n; ++i) {
		data[i] = f(time_step * i);
		freq_axis.push_back(i/options.time_range);
	}

	phys::stdVec_c transform = phys::FFT(data);

	phys::stdVec_d spectrum;
	for (auto& elem : transform) {
		elem /= n;
		spectrum.push_back(norm(elem));
	}

	phys::Viewer viewer;
	viewer.plot(freq_axis, spectrum);

	fprintf(stdout, "freq\t\tspectrum\n");
	for (int i = 0; i < n; ++i) {
		fprintf(stdout, "%f\t\t%f\n", freq_axis[i], spectrum[i]);
	}

	return 0;
}
