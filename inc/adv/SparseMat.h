#ifndef SPARSEMAT_HPP
#define SPARSEMAT_HPP

#include <adv/row.hpp>
#include <adv/mesh.hpp>
#include <staticMatrix.hpp>
#include <iomanip>

namespace phys{

	template<class T> class sparseMat;

	// transpose of a sparse matrix
	template<class T>
	const sparseMat<T> transpose(const sparseMat<T>&);

	// return the diagonal elements of a sparse matrix
	template<class T>
	const sparseMat<T> diagonal(const sparseMat<T>&);

	// Gauss-Seidel relaxation; x is overwriten with new iteration
	template<class T>
	void Gauss_Seidel(const sparseMat<T>& A, const std::vector<T>& f, std::vector<T>& x);

	// reverse Gauss-Seidel relaxation; x is overwriten with new iteration
	template<class T>
	void reverse_Gauss_Seidel(const sparseMat<T>& A, const std::vector<T>& f, std::vector<T>& x);

	// symmetric Gauss-Seidel relaxation; x is overwriten with new iteration
	template<class T>
	void symmetric_Gauss_Seidel(const sparseMat<T>& A, const std::vector<T>& f, std::vector<T>& x);

	template<class T>
	class sparseMat : public list<row<T>>{
	public:
		sparseMat(size_t n = 0);

		sparseMat(size_t n, const T& a);

		/**	Constructor using a mesh of triangles for the FEM */
		sparseMat(mesh<triangle>&);

		~sparseMat(){} // destructor

		const T operator() (size_t i, size_t j) const
		{
			return (*item[i])[j];
		}// read the (i,j)th element

		int row_number() const
		{
			return number;
		}// number of rows

		int column_number() const;

		int order() const
		{
			return std::max(row_number(), column_number());
		}// matrix order

		const sparseMat& operator+=(const sparseMat&);
		const sparseMat& operator-=(const sparseMat&);
		const sparseMat& operator*=(const T&);

		const sparseMat& truncate(double thresh);

		const sparseMat<T> factorize(double);
		const std::vector<T> forward_elimination(const std::vector<T>&)const;
		const std::vector<T> back_substitution(const std::vector<T>&)const;
		const std::vector<int> coarsen(double threshold = .05) const;
		const sparseMat<T> create_transfer();

		friend const sparseMat
			operator*<T>(const sparseMat&, const sparseMat&);

		friend const sparseMat transpose<T>(const sparseMat&);
		friend const sparseMat diagonal<T>(const sparseMat&);
		friend void Gauss_Seidel<T>(const sparseMat&, const std::vector<T>&, std::vector<T>&);
		friend void reverse_Gauss_Seidel<T>(const sparseMat&, const std::vector<T>&, std::vector<T>&);
		friend void symmetric_Gauss_Seidel<T>(const sparseMat&, const std::vector<T>&, std::vector<T>&);
	};

	// matrix plus matrix
	template<class T>
	const sparseMat<T> operator+(const sparseMat<T>&, const sparseMat<T>&);

	// matrix minus matrix
	template<class T>
	const sparseMat<T> operator-(const sparseMat<T>&, const sparseMat<T>&);

	// scalar times sparse matrix
	template<class T>
	const sparseMat<T> operator*(const T&, const sparseMat<T>&);

	// sparse matrix times scalar
	template<class T>
	const sparseMat<T> operator*(const sparseMat<T>&, const T&);

	// matrix times vector
	template<class T>
	const std::vector<T> operator*(const sparseMat<T>&, const std::vector<T>&);

	// sparse matrix times sparse matrix
	template<class T>
	const sparseMat<T> operator*(const sparseMat<T>&, const sparseMat<T>&);

	// vector, v, divided by main diagonal of a sparse matrix, A
	template<class T>
	const std::vector<T> operator/(	const std::vector<T>& v,
									const sparseMat<T>& A );

	// Jacobi relaxation; x is overwritten with new iteration.
	template<class T>
	void Jacobi(	const sparseMat<T>& A,
					const std::vector<T>& f,
					std::vector<T>& x );

	// ILU iterative lu solver.
	template<class T>
	void ILU(	const sparseMat<T>& A,
				const sparseMat<T>& L,
				const sparseMat<T>& U,
				const std::vector<T>& f,
				std::vector<T>& x	);

	// prints a sparse matrix to screen; includes zero elements.
	template<class T>
	std::ostream& operator<<(std::ostream&, const sparseMat<T>&);
}

#include <adv/SparseMat.inl>

#endif