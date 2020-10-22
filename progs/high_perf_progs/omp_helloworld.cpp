#include <iostream>

#include "libomp/omp.h"


int main () {

	//number of threads potentially available to OpenMP
	std::cout << "maximum num of threads: " << omp_get_max_threads() << std::endl;

	//number of logical processors available (may be different to the number of physical processors)
	std::cout << "num of logical processors: " << omp_get_num_procs() << std::endl;

	#pragma omp parallel
	{

		printf("printf: Hello from thread %d\n", omp_get_thread_num());

	}

	/*
	#pragma omp parallel
	{

		std::cout << "std::cout: Hello from thread " << omp_get_thread_num() << "\n";

	}
	*/

	return 0;
}
