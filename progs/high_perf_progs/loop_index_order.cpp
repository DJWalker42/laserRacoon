#include <iostream>

#include "Timer.h"

/*
 *  Contrived example to show the importance of understanding the
 *  underlying memory layout.
 *
 *
 *  Technical note: an array of arrays is created in contiguous memory
 *  i.e. an unbroken block.
 *  int a [rowCount][colCount]
 *
 *  Using new to create a two-dimensional array in a loop:
 *
 *  int** a = new int*[rowCount];
 *  for(int i = 0; i < rowCount; ++i)
 *  	a[i] = new int[colCount];
 *
 *  creates an array of pointers to arrays that are unlikely to be in a
 *  contiguous block of memory. Deleting this structure requires
 *  a corresponding loop and correct ordering; loop first to delete the
 *  nested arrays of a, then delete array a.
 *
 */


int main() {

	const int max_repeats {10};

	const int N {1000};
	const int N2 {N * N};

	//use new to ensure we create (large) memory on the heap
	double * A = new double [N2];
	double * B = new double [N2];
	double * C = new double [N2];

	phys::timer t;

	// i then j
	std::cout << "i then j\n";
	t.start();
	for (int repeat = 0; repeat < max_repeats; ++repeat) {
		for (int i = 0; i < N; ++i) {
			for (int j = 0; j < N; ++j) {
				C[j + N * i] = A[j + N * i] + B[j + N * i];
			}
		}
		C[0] -= 1.; //stops optimisations from removing the repeat loop
	}

	t.stop();

	t.display();

	//j then i
	std::cout << "j then i\n";
	t.start();
	for (int repeat = 0; repeat < max_repeats; ++repeat) {
		for (int j = 0; j < N; ++j) {
			for (int i = 0; i < N; ++i) {
				C[j + N * i] = A[j + N * i] + B[j + N * i];
			}
		}
		C[0] -= 1.; //stops optimisations from removing the repeat loop
	}
	t.stop();


	t.display();

	//not really needed as OS will clean up on exit but shows syntax
	delete[] A;
	delete[] B;
	delete[] C;


	return 0;
}
