#include <LinearSolvers.h>
#include <iostream>
using namespace phys;


int main(){

	size_t N = 3;

	mat A;

	stdVec_d row1(N), row2(N), row3(N);

	row1[0] = 25; row1[1] = 15; row1[2] = -5;
	row2[0] = 15; row2[1] = 18; row2[2] = 0;
	row3[0] = -5; row3[1] = 0 ; row3[2] = 11;

	A.append_row(row1);
	A.append_row(row2);
	A.append_row(row3);

	mat B;
	stdVec_d rowA(2), rowB(2), rowC(2);
	rowA[0] = 40; rowA[1] = -30;
	rowB[0] = 51; rowB[1] = -15;
	rowC[0] = 28; rowC[1] =  16;

	B.append_row(rowA); 
	B.append_row(rowB);
	B.append_row(rowC);

	LinearSolver* chol = new Cholesky( A );
	
	mat X = chol->solve(B);

	std::cout << "Cholesky\n";
	std::cout << X << "\n";

	row1[0] = 12; row1[1] = -51; row1[2] =   4;
	row2[0] =  6; row2[1] = 167; row2[2] = -68;
	row3[0] = -4; row3[1] = 24 ; row3[2] = -41;

	A.clear();

	A.append_row(row1);
	A.append_row(row2);
	A.append_row(row3);

	std::cout << "Orignal matrix:\n";
	std::cout << A << "\n";

	QR qr(A);

	std::cout << "QR decomp -> Q * R =\n";
	std::cout << qr.getQ() * qr.getR() << "\n";

	LU lu(A);

	std::cout << "LU decomp -> L * U =\n";
	std::cout << lu.getL() * lu.getU() << "\n";

	return 0;
}
