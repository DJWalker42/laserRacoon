#include "Fourier.h"
#include "Visualise.h"
#include "error_functions.h"

#include <cmath>

//doesn't exist on Windows machines
#include <unistd.h>  //getopt, optarg, optopt

struct options_t {
	int freq;   // test function: cos(2*pi*freq*t) * exp(-pi * t * t)
	int data_length; // number of "samples" for the test function
};

//an option letter followed by a colon means option requires a value
#define OPTION_FMT "f:n:h"
#define USAGE_FMT "%s [-f <frequency>] [-n <number of samples>\n]"

int freq; //global for convenience, avoids passing as argument to test function

double f(double x) {
	return cos(2. * M_PI * freq * x) * exp(-M_PI * x * x);
}

int main (int argc, char ** argv) {

	//default option values
	options_t options;
	options.freq = 3;
	options.data_length = 64;

	int c;
	while ((c = getopt(argc, argv, OPTION_FMT)) != -1) {
	    switch (c) {
	    case 'f':
	    	options.freq = atoi(optarg);
	    	break;
	    case 'n':
	    	options.data_length = atoi(optarg);
	    	break;
	    case '?':
	        if (optopt == 'f') {
	            phys::errExit("Option -f requires a value");
	        } else if (optopt == 'n') {
	        	phys::errExit("Option -n requires a value");
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

	freq = options.freq;
	int n {options.data_length};//for shorthand

	double t_min {-2.};
	double t_max {2.};
	double t_range = t_max - t_min;

	double step = t_range / (n);

	phys::stdVec_d time_axis(n);

	for(int i = 0; i < n; ++i) {
		time_axis[i] = t_min + i * step;
	}

	phys::stdVec_c data(n);

	for (int i = 0; i < n; ++i) {
		data[i] = f(time_axis[i]); //assigns to real part
	}

	phys::stdVec_c transform = phys::FFT(data);

	for (auto& elem : transform) {
		elem /= n;
		elem = conj(elem);
	}

	phys::stdVec_d spectrum;

	for (auto& elem : transform) {
		spectrum.push_back(norm(elem));
	}

	phys::stdVec_d freq_axis;
	for (int x = 0; x < n; ++x) {
		freq_axis.push_back(x / t_range);
	}

	phys::visual::Viewer viewer;
	viewer.withLines();
	viewer.plot(freq_axis, spectrum);

	return 0;
}
