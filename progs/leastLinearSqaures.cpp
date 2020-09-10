/*
 *  Program to solve the least linear squares fit for the example data
 *  	(1,2), (2,1), (3,3), (4,6)
 *  as explained in section 3.2.1 of the Computation Physics book.
 */

#include "LinearSolvers.h"


int main (int argc, char ** argv) {
	const size_t n = 2;
	phys::mat A;
	phys::stdVec_d row1(n) , row2(n);
	row1[0] = 4; row1[1] = 10;
	row2[0] = 10; row2[1]= 30;

	A.append_row(row1);
	A.append_row(row2);

	phys::Cholesky cholesky(A);

	phys::stdVec_d b {12,37};

	phys::stdVec_d x = cholesky.solve(b);

	std::cout << x << std::endl; //should output: -0.5 1.4

	return 0;
}
