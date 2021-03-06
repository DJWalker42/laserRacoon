#ifndef DJW_LINEARSOLVERS_HPP
#define DJW_LINEARSOLVERS_HPP

/*	The source code for the spline interpolation, including the Band Matrix definition,
 was read and modified from the following web address June 2014:
 http://kluge.in-chemnitz.de/opensource/spline/
 Author: Tino Kluge.
 */

#include "DynVector.h"
#include "DynMatrix.h"

#include <cstdint>

namespace phys {

/**	Virtual base class for matrix factorisation methods */
class LinearSolver {
protected:
	//Default constructor
	LinearSolver();
	//Constructor: A is the matrix you want to factor
	LinearSolver(const mat_d &A);
public:
	virtual ~LinearSolver() {
	}
private:
	/* No Copy */
	LinearSolver(const LinearSolver&);
	/* No Copy */
	const LinearSolver& operator=(const LinearSolver&);
public:
	//interface functions
	//decompose the matrix
	virtual void decompose(const mat_d &A = mat_d())=0;
	//solve a one dimensional system
	virtual stdVec_d solve(const stdVec_d &b) const = 0;
	//solve a multi-dimensional system
	virtual mat_d solve(const mat_d &B) const =0;
protected:
	mat_d m_factor;
};

// Cholesky Factorisation
class Cholesky: public LinearSolver {
public:
	Cholesky();
	Cholesky(const mat_d &A);
	void decompose(const mat_d &A = mat_d());
	/**
	 Solve a linear system A*x = b, using the previously computed
	 Cholesky factorisation of A: L*L'.

	 @param  b	A single (column) vector with as many elements as rows in A.
	 @return x	so that L*L'*x = b.
	 */
	stdVec_d solve(const stdVec_d &b) const;
	/**
	 Solve a linear system A*X = B, using the previously computed
	 Cholesky factorisation of A: L*L'.

	 @param  B	A Matrix with as many rows as A and any number of columns.
	 @return X	so that L*L'*X = B.
	 */
	mat_d solve(const mat_d &B) const;
	const mat_d getL() const {
		return m_factor;
	}
	bool is_spd() const {
		return m_spd;
	}
private:
	uint32_t m_dim;
	bool m_spd;
};

// QR Factorisation
class QR: public LinearSolver {
public:
	QR();
	QR(const mat_d &A);
	void decompose(const mat_d &A = mat_d());
	/**
	 Least squares solution of A*x = b
	 @param b    m-length vector
	 @returns	n-length vector that minimises the two norm of Q*R*X-B.
	 If B is non-conformant, or if QR.is_full_rank() is false,
	 the routine returns a null vector.
	 */
	stdVec_d solve(const stdVec_d &b) const;
	/**
	 Least squares solution of A*X = B
	 @param B    m x k Array (must conform).
	 @return X	n x k Array that minimises the two norm of Q*R*X-B.
	 If B is non-conformant, or if QR.is_full_rank() is false,
	 the routine returns a null matrix
	 */
	mat_d solve(const mat_d &B) const;

	//query functions

	/**
	 Generate and return the (economy-sized) orthogonal m_factor
	 @return     Q the (ecnomy-sized) orthogonal m_factor (Q*R=A).
	 */
	const mat_d getQ() const;
	/**
	 Return the upper triangular m_factor, R, of the QR factorisation
	 @return     R
	 */
	const mat_d getR() const;
	/**
	 Retreive the Householder vectors from the QR decomposition
	 @returns lower trapezoidal matrix whose columns define the reflections
	 */
	const mat_d getHouseholder() const;
	/**
	 Flag to denote if the matrix is of full rank or not.
	 @return true if matrix is full rank, false otherwise.
	 */
	bool is_full_rank() const;
private:
	uint32_t m_rows;
	uint32_t m_cols;
	stdVec_d m_Rdiag;
};

class LU: public LinearSolver {
public:
	LU();
	LU(const mat_d &A);
	/* Applies a version of the Doolittle algorithm; L_ii = 1.*/
	void decompose(const mat_d &A = mat_d());
	/**
	 Solve A*x = b, where x and b are vectors of length equal
	 to the number of rows in A.

	 @param  b   a vector of length equal to the first dimension of A.
	 @return x	a vector so that L*U*x = b(piv), if B is nonconformant,
	 returns a null vector.
	 */
	stdVec_d solve(const stdVec_d &b) const;
	/**
	 Solve A*X = B
	 @param  B   A Matrix with as many rows as A and any number of columns.
	 @return     X so that L*U*X = B(piv,:), if B is nonconformant, returns
	 a null matrix.
	 */
	mat_d solve(const mat_d &B) const;

	//query functions

	/**
	 Return lower triangular m_factor
	 @return     L
	 */
	mat_d getL() const;
	/**
	 Return upper triangular m_factor
	 @return     U portion of LU factorization.
	 */
	mat_d getU() const;
	/**
	 Return pivot permutation vector
	 @return     pivot
	 */
	stdVec_u getPivot() const;
	/**
	 Compute determinant using LU factors.
	 @return determinant of A, or 0 if A is not square.
	 */
	double det() const;
	/**
	 Is the matrix nonsingular?
	 @return     true if upper triangular m_factor U (and hence A)
	 is nonsingular, false otherwise.
	 */
	bool is_nonsingular() const;

private:
	// swaps the elements in b as per the pivot vector and returns the modified vector
	stdVec_d permute_copy(const stdVec_d &b) const;
	// swaps the rows in B as per the pivot vector and returns the modified matrix
	mat_d permute_copy(const mat_d &B) const;

private: //LU specific representation
	uint32_t m_rows;
	uint32_t m_cols;
	stdVec_u m_pivot;	//row i has been swapped with pivot[i].
	int m_pivsign;	//used in computation of the determinant
};

//	Band matrix solver - not a derived class of LinearSolver
//	implements its own specialised LU decomposition. NOTE: This does
//	not implement partial pivoting of rows - this is ok if the matrix
//	is diagonally dominant but may break if the matrix diagonal values
//  tend to zero (matrix -> singular).
class Band_matrix: public matrix<double> {
public:
	// Constructor
	Band_matrix(uint32_t dim, uint32_t nu, uint32_t nl);

	// public access operator
	//-- here i and j refer to the row and column indices of the entire matrix.
	// returns the corresponding element in the (reduced) band matrix.
	double operator ()(int i, int j) const;      // read
	double& operator ()(int i, int j);            // write

	// public interface functions.
	void lu_decompose();
	std::vector<double> r_solve(const std::vector<double> &b) const;
	std::vector<double> l_solve(const std::vector<double> &b) const;
	std::vector<double> lu_solve(const std::vector<double> &b);
private:
	//the addtional inverse diagonal is stored in last row vector in the band matrix.
	//we access the values using these read & write functions
	double saved_diag(int i) const;
	double& saved_diag(int i);
private:
	const int m_num_of_super;
	const int m_num_of_sub;
	const int m_dimensions;
	const uint32_t m_last_row_idx;
	bool m_is_lu_decomposed;
};
}

#endif
