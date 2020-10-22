#include "Timer.h"
#include "Visualise.h"
#include "Storage.h"

#include <numeric> //std::accumulate
#include <sstream>
#include <libomp/omp.h>
#include "phys_accumulate.h"

int main () {

	omp_set_dynamic(0); //disable dynamic "teams", ensures we get the number of threads we want

	//On Intel processors this number will be 2*number-of-cores, on AMD just the number-of-cores (?)
	int max_num_threads {omp_get_max_threads()};
	std::cout << "Max number of threads: " << max_num_threads << std::endl;
	std::cout << "Max number of processors: " << omp_get_num_procs() << std::endl;

	const int max_repeats {10};

	// processor on my machine Intel quad-core 2.2GHz, used to compute a "model" time for one thread (serial)
	const double processor_freq {2.2e9};//you should change this value for your machine

	/*
	 *  There's a bug (feature?) in Storage class whereby we have to set the names
	 *  for each dependent variable we are going to store else we get a runtime error.
	 *  Here we need to match the loop governing the number of threads to use when
	 *  calling the OMP version of accumulate.
	 */
	std::vector<std::string> names;
	names.push_back("model");
	names.push_back("std");
	for (int threads = 1; threads <= max_num_threads; threads *= 2) {
		std::stringstream ss;
		ss << "omp_T" << threads;
		names.push_back(ss.str());
	}
	names.push_back("phys");

	phys::Storage<double> container("n", names); //"n" is the independent (x-axis) name

	phys::timer outer;
	outer.start();

	for (int n = 1000000; n <= 50000000; n += 1000000) {
		std::vector<double> v(n, 1.); //vector size n all elements 1. => sum == n
		phys::timer t;

		std::vector<double> times; //temporary store of the timings

		//compute the "model" T1 (serial) timing for accumulate, based on 2 clock-cyles per element
		times.push_back(2 * n / processor_freq);

		//serial std::accumulate
		t.start();
		for (int repeat = 0; repeat < max_repeats; ++repeat) {
			auto sum = std::accumulate(std::begin(v), std::end(v), 0. );

			/*
			 *  Compiler optimisations may do some interesting things here if
			 *  we don't do anything with the result of std::accumulate.
			 *  You can test this by commenting out the output line below,
			 *  compiling with -O2, ignore the warning about 'sum', and looking
			 *  at the timing results.
			 *
			 *
			 *  (the if condition and output add a constant amount of
			 *  time to the overall timing of the repeats but this becomes
			 *  decreasingly significant as n increases)
			 */

			if (repeat == 0) std::cout << "std sum =\t" << sum << "\n";
		}
		t.stop();
		times.push_back(t.get() / max_repeats);

		//parallel OMP accumulate, num threads 1,2,4,8 ..., max_num_threads
		for(int threads = 1; threads <= max_num_threads; threads *= 2) {
			omp_set_num_threads(threads);

			t.start();
			for (int repeat = 0; repeat < max_repeats; ++repeat) {
				auto sum = phys::omp_accumulate(v, 0.);
				if (repeat == 0) std::cout << "omp T" << threads << ": sum =\t" << sum << "\n";
			}
			t.stop();

			times.push_back(t.get() / max_repeats);
		}


		t.start();
		for (int repeat = 0; repeat < max_repeats; ++repeat) {
			auto sum = phys::accumulate_p(std::begin(v), std::end(v), 0.);
			if (repeat == 0) std::cout << "phys sum =\t" << sum << "\n";
		}
		t.stop();
		times.push_back(t.get() / max_repeats);

		container.store(double(n), times);
	}
	outer.stop();
	outer.display();

	container.write("./accumulate.log", true); //write headers

	phys::Viewer viewer(1600, 900); //width, height in pixels
	viewer.plot(container);

	return 0;
}
