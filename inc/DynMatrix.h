#ifndef PHYSMATRIX_HPP
#define PHYSMATRIX_HPP

#include "DynVector.h"
#include "Maths.h"

namespace phys{

	enum dimension {HEIGHT, WIDTH}; 

	//------------------------------------------------------------------------
	//	Matrix class
	//------------------------------------------------------------------------
	//  mat[i] ALWAYS gives ith row vector with the understanding that
	//	the jth column vector of A is equivalent to the jth row vector of A^T
	//
	// Note: The trans flag merely indicates we wish to treat the matrix as its 
	// transpose not actually physically transpose the matrix. The non-member
	// function transpose does this.
	template<class T>
	class matrix {			
	public:
		//default constructor
		matrix( size_t m = 0 );

		// constructor for m-by-n matrix of zeros	
		matrix(	size_t m, size_t n, bool is_T = false);

		// constructor for m-by-n matrix of val
		matrix( size_t m, size_t n, const T& val, bool is_T = false);

		// constructor for m-by-row.size() matrix of repeated row
		matrix( size_t m, const std::vector<T>& row, bool is_T = false);

		// constructor using a vector of vectors.
		matrix( const std::vector<std::vector<T>>& vv, bool is_T = false);

		/*	constructor using a single vector with the number of rows
		and columns specified (here v.size() must equal rows * cols) */
		matrix( const std::vector<T>& v, size_t rows, size_t cols, bool is_T = false);

		// Copy constructor
		matrix( const matrix& mat );

		matrix& operator=( const matrix& ); 

		/*	Element read access operator: (row, col) */
		const T& operator()(size_t i, size_t j) const;

		/*	Element write access operator: (row, col) */
		T& operator()(size_t i, size_t j); 

		/* Read access to the ith row/column vector of the matrix - depends on trans flag */ 
		const std::vector<T>& operator[](size_t i) const;

		/* Write access to the ith row/column vector of the matrix - depends on trans flag */
		std::vector<T>& operator[](size_t i); 

		/* Select a sub-matrix of current matrix - this is read only*/
		const matrix operator()( size_t imin, size_t imax, size_t jmin, size_t jmax ) const;

		/* If the value assigned is 1 creates the identity matrix, else broadcasts the value to all elements */
		const matrix& operator=(const T&);

		const matrix& operator+=(const matrix&);
		const matrix& operator-=(const matrix&);

		const matrix& operator+=(const T&);
		const matrix& operator-=(const T&);
		const matrix& operator*=(const T&);
		const matrix& operator/=(const T&);

		void append_row(const std::vector<T>&);
		void append_col(const std::vector<T>&);

		/* Calls std::vector clear() on the data array*/
		void clear();
		/* Calls std::vector empty() on the data array*/
		bool empty();

		const matrix& swap_rows(size_t, size_t, bool = false);
		const matrix& swap_cols(size_t, size_t, bool = false);

		const matrix& resize(size_t, const std::vector<T>& = std::vector<T>(), dimension which = HEIGHT);
		/*	Resizes matrix: if reducing the dimensions thats all it does, if increasing the dimensions the
			matrix is padded with zeros (default) or the value you supply */
		const matrix& pad_matrix(size_t, size_t, const T& = T());

		size_t size() const { return data_array.size(); }

		bool empty() const { return data_array.empty(); }

		//Query the transpose state
		bool is_trans() const {return trans;}

		//Switch the transpose state
		void transpose(){ trans = (trans) ? false : true; }

	private:
		std::vector<std::vector<T>> data_array; //!< The container for the matrix elements
		bool trans;								//!< is the matrix its transpose or not?
	};

	//***************************************
	//matrix non-member operator functions
	//***************************************

	template<typename T>
	const matrix<T> operator+(const matrix<T>&m1, const matrix<T>&m2);

	template<typename T>
	const matrix<T> operator-(const matrix<T>&m1, const matrix<T>&m2);

	template<typename T>
	const matrix<T> operator+(const T&a, const matrix<T>&m);

	template<typename T>
	const matrix<T> operator+(const matrix<T>&m, const T&a);

	template<typename T>
	const matrix<T> operator-(const T&a, const matrix<T>&m);

	template<typename T>
	const matrix<T> operator-(const matrix<T>&m, const T&a);

	template<typename T>
	const matrix<T> operator+(const matrix<T>&m); 

	template<typename T>
	const matrix<T> operator-(const matrix<T>&m);

	//multiplication scalar by matrix
	template<class T>
	const matrix<T> operator*(const T&a, const matrix<T>&m);

	//multiplication matrix by scalar
	template<class T>
	const matrix<T> operator*(const matrix<T>&m, const T&a);

	//division matrix by scalar
	template<class T>
	const matrix<T> operator/(const matrix<T>&m, const T&a);

	//vector times matrix (row vector v by column vectors of m)
	template<class T>
	const std::vector<T> operator*(const std::vector<T>&v, const matrix<T>&m);

	//matrix times vector (row vectors of m by column vector v)
	template<class T>
	const std::vector<T> operator*(const matrix<T>&m, const std::vector<T>&v);

	/*	Note the following about matrix multiplication; ' means transpose:
		C(M,K)  = A(M,N) B(N,K)
		C'(K,M) = B'(K,N) A'(N,M);
	*/

	// matrix multiplication
	template<class T>
	const matrix<T> operator*(const matrix<T>&m1, const matrix<T>&m2);

	/*	For performance purposes, this function assumes you have correctly set up
		the second matrix. That is if you want the result C = A * B, then B 
		has to be transposed e.g. B' computed from transpose(B); the flag set by B.transpose()
		will not work if using this function directly. 
		If you want C = A * B' then pass B without transposing. A matrix, B, multiplied by its transpose
		is simply computed by calling matmul(B,B)	*/
	template<class T>
	const matrix<T> matmul(const matrix<T>& m1, const matrix<T>& m2);

	/*	As matmul that returns a matrix but here overwrites m1 with the result.
		Saves on memory. */
	template<class T>
	void _matmul(matrix<T>& m1, const matrix<T>& m2);

	/** Transpose the elements of a matrix; returns the transposed matrix.
		Note that this funtion ignores the trans flag of the input matrix,  m. 
		Resultant matrix is not copy constructed+ from m and thus has a trans flag
		of false, regardless of m's trans flag. +(except when the matrix is empty
		or a single vector; in either case the flag is redundant anyway).	*/
	template<class T>
	const matrix<T> transpose(const matrix<T>&m);

	//transpose a matrix "in-place" saves on memory - if the matrix is not square could be computationally expensive.
	template<class T>
	void _transpose( matrix<T>& m );

	//creates a dim-by-dim identity matrix
	//we have to give the complete type specification in the call. 
	template<typename T>
	const matrix<T> identity(size_t dim);

	//output operator for matrices.
	template<class T>
	std::ostream & operator<<(std::ostream& os, const matrix<T>& m);

	//swap row i with row j im matrix m
	template<class T>
	const matrix<T> swap_rows(const matrix<T>&m, size_t i, size_t j);

	//swap column i with column j in matrix m 
	template<class T>
	const matrix<T> swap_cols(const matrix<T>&m, size_t i, size_t j);

	//resize the height of a matrix 
	template<class T>
	const matrix<T> resize(const matrix<T>&m, size_t resize_size, const std::vector<T>& val = std::vector<T>(),
		dimension which = HEIGHT);

	//resize width and height of a matrix
	template<class T>
	const matrix<T> pad_matrix(const matrix<T>&m, size_t rows, size_t cols, const T& val = T());

	// read a matrix from a text file; To call we have to give the complete type specfication
	// e.g. mat A; A = read_matrix<double>(fn); otherwise complier cannot deduce the type.
	template<class T>
	matrix<T> read_matrix( const std::string& filename, bool transpose = false );

	typedef matrix<double>	mat;
	typedef matrix<int>		mat_i;
	typedef matrix<size_t>	mat_s;
}

#include "DynMatrix.inl"

#endif
