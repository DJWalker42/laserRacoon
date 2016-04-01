#ifndef STATICMATRIX_HPP
#define STATICMATRIX_HPP

#include "StaticVector.h"

namespace phys{

	//Forward declaration of the staticMat class
	//- used in declaration of the forward declaration output stream operator
	template<class T, size_t M, size_t N > class staticMat;

	template<class T, size_t M, size_t N>
	std::ostream& operator<<(std::ostream& os, const staticMat<T,M,N>& mat);

	/** Static Matrix class
		Convention here is M rows by N columns
	*/
	template<class T, size_t M, size_t N>
	class staticMat : public staticVec<staticVec<T,N>, M >{
	public:
		/********** Constructors *************/

		/* Default constructor - creates an M x N matrix of zeros */
		staticMat();

		/* Constructor using a single vector - vector is repeated in all M rows */
		staticMat(	const staticVec<T,N>& v, 
					bool t = false	);

		/* Special constructor for a staticMat with two rows (M) (v,u) */
		staticMat(	const staticVec<T,N>& v, 
					const staticVec<T,N>& u, 
					bool t = false	);

		/* Special constructor for a staticMat with three rows (v,u,w) */
		staticMat(	const staticVec<T,N>& v, 
					const staticVec<T,N>& u, 
					const staticVec<T,N>& w, 
					bool t = false	);

		/*	Special constructor with T argument; only diagonal elements get the value;
			typically used to create the identity matrix (i.e. val = 1) */
		staticMat(	const T& val,
					bool t = false	);						

		/*	Constructor from base class object */
		staticMat(	const staticVec<staticVec<T,N>, M>& vv,
					bool t = false );

		/* Copy constructor */
		staticMat(	const staticMat&  );

		/********** Operators *************/ 

		/* Copy assignment operator */
		const staticMat& operator=( const staticMat& );

		/*	FORTRAN like access operator -- read
			Alternative to [i][j] */
		const T& operator()(size_t i, size_t j) const;

		/*	FORTRAN like access operator -- write 
			Alternative to [i][j] */
		T& operator()(size_t i, size_t j); 

		//special assignment to set diagonal elements to val; typically used to create the identity matrix
		const staticMat& operator=(const T& val);

		//current matrix (lhs) times matrix rhs in that order, result stored to current matrix.
		const staticMat& operator*=(const staticMat& rhs);

		//current matrix times scalar
		const staticMat& operator*=(const T&);

		//current matrix divided by scalar
		const staticMat& operator/=(const T&);

		/*	In client code ignore the bool argument and call these as normal e.g. m.l_one_norm().
			If trans flag set true, i.e. we are treating m as its transpose, it will automatically
			switch to the l_inf_norm() function or vice versa. */

		//l one norm of the current matrix. Ignore bool argument and call as m.l_one_norm() 
		const T l_one_norm(bool = false) const;

		//l infinity norm of the current matrix. Ignore bool argument and call as m.l_inf_norm()
		const T l_inf_norm(bool = false) const;

		//query the transpose flag
		bool is_trans() const {return m_trans;}

		/* Switches transpose state */
		void transpose(){ m_trans = (m_trans) ? false : true;}

		/*	Output */
		friend std::ostream& operator<< <T, M, N>(std::ostream&, const staticMat&);
	private:
		bool m_trans;		//!< is the matrix to be treated as it's transpose or not?
	};

	/* Non-member operators and functions */

	//staticMat times scalar
	template<class T, size_t M, size_t N>
	const staticMat<T,M,N> operator*(const staticMat<T,M,N>& m, const T& a);

	//scalar times staticMat
	template<class T, size_t M, size_t N>
	const staticMat<T,M,N> operator*(const T& a, const staticMat<T,M,N>& m);

	//staticMat divded by scalar
	template<class T, size_t M, size_t N>
	const staticMat<T,M,N> operator/(const staticMat<T,M,N>& m, const T& a);

	// row vector times matrix.
	template<class T, size_t M, size_t N>
	const staticVec<T,N> operator*(	const staticVec<T,M>& v, 
									const staticMat<T,M,N>& m);

	// matrix times column vector
	template<class T, size_t M, size_t N>
	const staticVec<T,M> operator*(	const staticMat<T,M,N>& m, 
									const staticVec<T,N>& v	);

	// matrix multiplication
	template<class T, size_t M, size_t N, size_t K>
	const staticMat<T,M,K> operator*(	const staticMat<T,M,N>&m1, 
										const staticMat<T,N,K>&m2	);

	/*	If using this function directly note that it assumes you have correctly set up
		the second matrix. That is if you want the result C = A * B, then B has to be 
		transposed i.e. B' computed from transpose(B); this is done for performance 
		reasons. Note that the flag set by B.transpose() will not work. If you want 
		C = A * B' then pass B without transposing. A matrix mulitplied by its 
		transpose is then simply computed by calling matmul(B, B). */
	template<class T, size_t M, size_t N, size_t K>
	const staticMat<T,M,K> matmul(	const staticMat<T,M,N>& m1, 
									const staticMat<T,K,N>& m2	); //note dimension order in m2

	/*	As matmul that returns a matrix but here overwrites m1 with the result. Saves on memory.*/
	template<class T, size_t M, size_t N, size_t K>
	void _matmul(staticMat<T,M,N>& m1, const staticMat<T,N,K>& m2);

	//transpose the elements of a matrix; returns the transposed matrix.
	template<class T, size_t M, size_t N>
	const staticMat<T,N,M> transpose(const staticMat<T,M,N>&m);

	// in place transpose; saves on memory; square matrices only.
	template<class T, size_t M, size_t N>
	void _transpose( staticMat<T,M,N>& m );

	// estimates upper bound for l_2 norm; l_one and l_inf norm always postive.
	template<class T, size_t M, size_t N>
	const T l_2_norm(const staticMat<T,M,N>& m);

	//output stream for static matrices
	template<class T, size_t M, size_t N>
	std::ostream& operator<<(std::ostream& os, const staticMat<T,M,N>& m);

	//type defintion for 2 x 2 matrices
	typedef staticMat<double, 2, 2> staticMat2;
	//type defintion for 3 x 3 matrices
	typedef staticMat<double, 3, 3> staticMat3;

	//Determinant of a 2x2 matrix
	double det(const staticMat2& A);

	// Inverse of a 2x2 matrix by Kramer's formula.
	const staticMat2 inverse(const staticMat2& A);
}

#include "StaticMatrix.inl"

#endif
