/*
 *  Program to solve the least linear squares fit for the example data
 *  	(1,2), (2,1), (3,3), (4,6)
 *  as explained in section 3.2.1 of the Computation Physics book.
 *
 *  Hint for exercise 3:
 *
 *  n == number of data pairs
 *  m == order of fitting polynomial
 *
 *  S = sum i=0 -> n [u^2] => dS/du = sum i=0 -> n [2u]
 *  u = yi - a0 - a1.xi - ... - am.xi^m
 *  => du/da0 = -1, du/da1 = -xi, ..., du/dam = -xi^m  (partial derivatives)
 *
 *  dS/dak = dS/du . du/dak (chain rule)
 *
 *  RHS vector 'b' then consists of terms with no ak product, elements
 *  increasing in order up to m.
 *
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
