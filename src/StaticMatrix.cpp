#include "StaticMatrix.h"

namespace phys {

double det(const staticMat2 &A) {
	return A(0, 0) * A(1, 1) - A(0, 1) * A(1, 0);
} //determinant of a 2x2 matrix

const staticMat2 inverse(const staticMat2 &A) {
	point2 row0 =
			A.is_trans() ?
					point2(A(1, 1), -A(1, 0)) : point2(A(1, 1), -A(0, 1));
	point2 row1 =
			A.is_trans() ?
					point2(-A(0, 1), A(0, 0)) : point2(-A(1, 0), A(0, 0));
	return staticMat2(row0, row1) / det(A);
} // inverse of a 2x2 matrix by Kramer's formula.

} //namespace
