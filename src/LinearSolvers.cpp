#include <cassert>
#include <algorithm>
#include <iostream>

#include "LinearSolvers.h"
#include "Maths.h"

/*	The source code for the spline interpolation, including the Band Matrix definition,
	was read and modified from the following web address June 2014:			
		http://kluge.in-chemnitz.de/opensource/spline/ 
	Author: Tino Kluge.

	The matrix factorisation algorithms have been modified from the JAMA headers.
*/

namespace phys{
	// -------------------------------------------------------
	// Matrix Factorization Implementation
	// -------------------------------------------------------

	LinearSolver::LinearSolver() : m_factor()
	{}

	LinearSolver::LinearSolver(const mat& A) : m_factor(A) 
	{}

	// -------------------------------------------------------
	// Cholesky Implementation
	// -------------------------------------------------------

	Cholesky::Cholesky() : LinearSolver(), m_dim(0), m_spd(false)
	{}

	Cholesky::Cholesky(	const mat&A ) : LinearSolver(A),
						m_dim(A.size()),
						m_spd(true)
	{
		decompose();
	}

	void Cholesky::decompose( const mat&A )
	{
		if(!A.empty())
		{
			m_factor = A;
			m_dim = m_factor.size();
			m_spd = true;
		}
		//else m_factor created from A in constructor
		if(m_factor.empty()) return;

		//check A is square
		if( m_dim != m_factor[0].size())
		{
			std::cout << "Matrix not sqaure" << std::endl;
			m_spd = false;
			return;
		}
		//check A is symmetric
		for( size_t i = 0; i < m_dim; ++i )
		{
			for( size_t j = 0; j < m_dim; ++ j )
			{
				if( m_factor[i][j] != m_factor[j][i] )
				{
					std::cout << "Matrix not symmetric at " 
						<< i << ", " << j << std::endl;
					m_spd = false;
					return;
				}
			}
		}
		//create a local matrix to store the lower triangle m_factor
		//zero assigned to all elements by default.
		mat L(m_dim, m_dim);

		// Main loop.
		for (size_t j = 0; j < m_dim; j++) 
		{
			double d = 0.0;
			for (size_t k = 0; k < j; k++) 
			{
				double s = 0.0;
				for (size_t i = 0; i < k; i++) 
				{
					s += L[k][i]*L[j][i];
				}
				L[j][k] = s = (m_factor[j][k] - s)/L[k][k];
				d = d + s*s;
			}
			d = m_factor[j][j] - d;
			if( d < 0.0 ) //check A is positive definite
			{
				std::cout << "Matrix not positivie definite: fails at:" 
					<< j << ", " << j << "\n";
				m_spd = false;
				m_factor = L; //assign what we have computed to m_factor
				return;
			}
			L[j][j] = sqrt(d);
		}
		//assign L to m_factor.
		m_factor = L; 
		return;
	}

	stdVec_d Cholesky::solve(const stdVec_d &b) const
	{
		if(m_factor.empty() || !m_spd)
		{
			std::cout << "Matrix A not factorised; returning null\n";
			return stdVec_d();
		}
		//size_t m_cols = m_factor.size();
		if (b.size() != m_dim)
		{
			std::cout << "b is non-conformant; returning null\n";
			return stdVec_d();
		}
		stdVec_d x(b);
		// Solve L*y = b;
		for (size_t k = 0; k < m_dim; ++k) 
		{
			for (size_t i = 0; i < k; ++i) 
				x[k] -= x[i]*m_factor[k][i];
			x[k] /= m_factor[k][k];
	
		}
		// Solve L'*x = y;
		for (int k = m_dim - 1; k >= 0; --k) //using an int due to >= 0 condition.
		{
			for (size_t i = k + 1; i < m_dim; i++){
				x[k] -= x[i] * m_factor[i][k];
			}
			x[k] /= m_factor[k][k];
		}
		return x;
	}

	mat Cholesky::solve( const mat &B ) const
	{
		if(m_factor.empty() || !m_spd)
		{
			std::cout << "Matrix A factorisation failed; returning null\n";
			return mat();
		}

		mat X(B);
		if(B.is_trans()) _transpose(X); //in-place transpose (matrix always square)
		if (X.size() != m_dim )
		{
			std::cout << "B is non-conformant; returning null\n";
			return mat();
		}

		size_t nx = X[0].size(); //num of cols

		// Solve L*Y = B;
		for (size_t j=0; j< nx; ++j)
		{
			for (size_t k = 0; k < m_dim; ++k) 
			{
				for (size_t i = 0; i < k; ++i){ 
					X[k][j] -= X[i][j] * m_factor[k][i];
				}
				X[k][j] /= m_factor[k][k];
			}
		}

		// Solve L'*X = Y;
		for (size_t j = 0; j < nx; ++j)
		{
			for (int k = m_dim - 1; k >= 0; --k) //using an int due to >= 0 condition.
			{
				for (size_t i = k+1; i < m_dim; ++i){ 
					X[k][j] -= X[i][j] * m_factor[i][k];
				}
				X[k][j] /= m_factor[k][k];
			}
		}
		return X;
	}

	// ---------------------------------------------------------
	// QR Implementation
	// ---------------------------------------------------------

	QR::QR() : LinearSolver(), m_rows(0), m_cols(0), m_Rdiag(m_cols)
	{}

	QR::QR(	const mat&A	) : LinearSolver(A),
			m_rows(A.size()),
			m_cols(A[0].size()),
			m_Rdiag(m_cols)
	{
		decompose();
	}

	void QR::decompose (const mat &A)		
	{
		if(!A.empty()){
			m_factor = A; //assign argument A
			m_rows = m_factor.size();
			m_cols = m_factor[0].size();
			m_Rdiag.resize(m_cols); //values overwritten during algorithm
		}
		//else m_factor copied from A in constructor

		if(m_factor.empty()) return; //chk m_factor has been set

		// Main loop.
		for (size_t k = 0; k < m_cols; ++k) {
			// Compute 2-norm of k-th column without under/overflow.
			double nrm = 0;
			for (size_t i = k; i < m_rows; ++i) 
			{
				nrm = phys::maths::mag2D(nrm,m_factor[i][k]);
			}

			if (nrm != 0.0) 
			{
				// Form k-th Householder vector.
				if (m_factor[k][k] < 0) nrm = -nrm;
				for (size_t i = k; i < m_rows; i++) 
					m_factor[i][k] /= nrm;

				m_factor[k][k] += 1.0;

				// Apply transformation to remaining columns.
				for (size_t j = k + 1; j < m_cols; ++j) 
				{
					double s = 0.0; 
					for (size_t i = k; i < m_rows; ++i) 
						s += m_factor[i][k] * m_factor[i][j];

					s = -s / m_factor[k][k];

					for (size_t i = k; i < m_rows; ++i) 
						m_factor[i][j] += s * m_factor[i][k];						
				}
			}
			m_Rdiag[k] = -nrm;
		}
	}

	bool QR::is_full_rank() const		
	{
		for (size_t j = 0; j < m_cols; ++j) 
		{
			if (m_Rdiag[j] == 0)
			return false;
		}
		return true;
	}

	const mat QR::getHouseholder() const
	{
		mat H(m_rows, m_cols); //initialises zero everywhere

		for (size_t i = 0; i < m_rows; ++i) 
		{
			for (size_t j = 0; j < m_cols; ++j) 
			{
				if (i >= j) 
					H[i][j] = m_factor[i][j]; 
				//else { H[i][j] = 0.0;} upper triangle
			}
		}
		return H;
	}	

	const mat QR::getR() const
	{
		mat R(m_cols,m_cols); //intitialised with zeros everywhere
		for (size_t i = 0; i < m_cols; ++i) 
		{
			for (size_t j = 0; j < m_cols; ++j) 
			{
				if (i < j) 
					R[i][j] = m_factor[i][j];
				else if (i == j) 
					R[i][j] = m_Rdiag[i];
				//else { R[i][j] = 0.0; } -- lower triangle
			}
		}
		return R;
	}

	const mat QR::getQ() const
	{
		mat Q(m_rows,m_cols);
		for (int k = m_cols - 1; k >= 0; --k) //int due to >=0 condition
		{ 
			Q[k][k] = 1.0;
			for (size_t j = k; j < m_cols; ++j) 
			{
				if (m_factor[k][k] != 0) 
				{
					double s = 0.0;
					for (size_t i = k; i < m_rows; ++i)
						s += m_factor[i][k] * Q[i][j];
				
					s = -s / m_factor[k][k];

					for (size_t i = k; i < m_rows; ++i) 
						Q[i][j] += s * m_factor[i][k];						
				}
			}
		}
		return Q;
	}

	stdVec_d QR::solve(const stdVec_d &b) const
	{
		if(b.size() != m_rows) return stdVec_d();
		if( !is_full_rank() ) return stdVec_d();
		stdVec_d x(b);

		// Compute Y = transpose(Q)*b
		for (size_t k = 0; k < m_cols; ++k) 
		{
			double s = 0.0; 
			for (size_t i = k; i < m_rows; ++i) 
			{
				s += m_factor[i][k] * x[i];
			}
			s = -s / m_factor[k][k];
			for (size_t i = k; i < m_rows; ++i) 
			{
				x[i] += s * m_factor[i][k];
			}
		}
		// Solve R*X = Y;
		for (int k = m_cols - 1; k >= 0; --k) //int due to >= 0 condition
		{
			x[k] /= m_Rdiag[k];
			for (size_t i = 0; i < size_t(k); ++i) 
				x[i] -= x[k] * m_factor[i][k];				
		}

		/* return m_cols portion of X */
		stdVec_d x_(m_cols);
		for (size_t i = 0; i < m_cols; ++i)
			x_[i] = x[i];

		return x_;
	}

	mat QR::solve(const mat &B) const
	{
		if( B.size() != m_rows ) return mat(); 
		if( !is_full_rank() ) return mat();

		size_t nx = B[0].size();
		mat X(B);

		// Compute Y = transpose(Q)*B
		for (size_t k = 0; k < m_cols; ++k) 
		{
			for (size_t j = 0; j < nx; ++j) 
			{
				double s = 0.0; 
				for (size_t i = k; i < m_rows; ++i) 
					s += m_factor[i][k] * X[i][j];
				
				s = -s / m_factor[k][k];
				for (size_t i = k; i < m_rows; ++i) 
					X[i][j] += s * m_factor[i][k];
				
			}
		}
		// Solve R*X = Y;
		for (int k = m_cols - 1; k >= 0; --k) //int due to >= 0 condition
		{ 
			for (size_t j = 0; j < nx; ++j)
				X[k][j] /= m_Rdiag[k];
			
			for (size_t i = 0; i < size_t(k); ++i) 
				for (size_t j = 0; j < nx; ++j) 
					X[i][j] -= X[k][j] * m_factor[i][k];				
		}

		/* return m_cols x nx portion of X */
		mat X_(m_cols,nx);
		for (size_t i = 0; i < m_cols; ++i)
			for (size_t j = 0; j < nx; ++j)
				X_[i][j] = X[i][j];

		return X_;
	}

	// --------------------------------------------------------
	// LU Implementation
	// --------------------------------------------------------

	LU::LU() : LinearSolver(), m_rows(0), m_cols(0), m_pivot(m_rows)
	{}

	LU::LU(	const mat&A) : LinearSolver(A),
			m_rows(A.size()),
			m_cols(A[0].size()),
			m_pivot(m_rows)
	{
		decompose();
	}

	void LU::decompose( const mat& A )
	{
		if(!A.empty()){ //i.e. calling function with an argument
			m_factor = A;
			m_rows = m_factor.size();
			m_cols = m_factor[0].size();
			m_pivot.resize(m_rows); //values overwritten at start of algorithm.
		}
		//else m_factor copied from A in constructor
		if(m_factor.empty()) return;

		// Use a "left-looking", dot-product, Doolittle algorithm.
		for (size_t i = 0; i < m_rows; i++) 
			m_pivot[i] = i;
		
		m_pivsign = 1;
		stdVec_d LUcolj(m_rows);

		// Outer loop.
		for (size_t j = 0; j < m_cols; j++) 
		{
			// Take a copy of the j-th column.
			for (size_t i = 0; i < m_rows; i++) 
				LUcolj[i] = m_factor[i][j];
			
			// Apply previous transformations using dot product
			for (size_t i = 0; i < m_rows; i++) 
			{
				size_t kmax = std::min(i,j);
				double s = 0.0;
				for (size_t k = 0; k < kmax; k++) 
					s += m_factor[i][k] * LUcolj[k];
				
				m_factor[i][j] = LUcolj[i] -= s; 
			}

			// Find m_pivot and exchange if necessary.
			size_t p = j;
			for (size_t i = j + 1; i < m_rows; i++) 
			{
				if (fabs(LUcolj[i]) > fabs(LUcolj[p])) 
					p = i;
				
			}
			if (p != j) 
			{
				for (size_t k = 0; k < m_cols; k++) 
				{
					double t = m_factor[p][k]; 
					m_factor[p][k] = m_factor[j][k]; 
					m_factor[j][k] = t;
				}
				size_t kp = m_pivot[p]; 
				m_pivot[p] = m_pivot[j]; 
				m_pivot[j] = kp;
				m_pivsign = -m_pivsign;
			}

			// Compute multipliers.        
			if ((j < m_rows) && (m_factor[j][j] != 0.0)) 
			{
				for (size_t i = j + 1; i < m_rows; i++) 
					m_factor[i][j] /= m_factor[j][j];
			}
		}
	}

	bool LU::is_nonsingular() const {
		size_t jmax = std::min(m_cols,m_rows);
		for (size_t j = 0; j < jmax; j++) 
		{
			if (m_factor[j][j] == 0)
				return false;
		}
		return true;
	}

	mat LU::getL() const {
		mat L_;
		if( m_rows > m_cols)	L_ = mat(m_rows,m_cols);
		else		L_ = mat(m_rows,m_rows);//m_rows <= m_cols

		for (size_t i = 0; i < m_rows; ++i) 
		{
			for (size_t j = 0; j < std::min(i,m_cols); ++j) 
				L_[i][j] = m_factor[i][j];

			L_[i][i] = 1.0;
		}		
		return L_;
	}

	mat LU::getU() const {
		mat U_;
		if (m_rows > m_cols)	U_ = mat(m_cols,m_cols);
		else		U_ = mat(m_rows,m_cols);//m_rows <= m_cols
		size_t imax = std::min(m_rows,m_cols);
		for (size_t i = 0; i < imax; ++i) 
			for (size_t j = i; j < m_cols; ++j) 
				U_[i][j] = m_factor[i][j];

		return U_;
	}

	stdVec_s LU::getPivot() const 
	{
		return m_pivot;
	}

	double LU::det() const {
		if (m_rows != m_cols) return double(0);
		double d = double(m_pivsign);

		for (size_t j = 0; j < m_cols; j++) 
			d *= m_factor[j][j];

		return d;
	}

	stdVec_d LU::solve(const stdVec_d &b) const
	{
		/* Dimensions: A(m_rows,m_cols), x(m_cols), b(m_rows) */
		if (b.size() != m_rows) return stdVec_d();
		if (!is_nonsingular()) return stdVec_d();

		stdVec_d x = permute_copy(b);

		// Solve L*Y = B(piv)
		for (size_t k = 0; k < m_cols; ++k)
			for (size_t i = k + 1; i < m_cols; ++i)
				x[i] -= x[k] * m_factor[i][k];

		// Solve U*X = Y;
		for (int k = m_cols - 1; k >= 0; --k) //int due to >= 0 condition
		{
			x[k] /= m_factor[k][k];
			for (size_t i = 0; i < size_t(k); i++)
				x[i] -= x[k] * m_factor[i][k];
		}
		return x;
	}

	mat LU::solve (const mat &B) const
	{
		/* Dimensions: A(m_rows,m_cols), X(m_cols,k), B(m_rows,k) */   
		if (B.size() != m_rows) {std::cout << "LU: B non-conformant \n";return mat();}
		if (!is_nonsingular()) {std::cout << "LU: A non-singular \n"; return mat();}

		// Copy right hand side with pivoting applied
		size_t nx = B[0].size();
		mat X = permute_copy(B);

		// Solve L*Y = B(piv,:)
		for (size_t k = 0; k < m_cols; ++k) 
			for (size_t i = k + 1; i < m_cols; ++i) 
				for (size_t j = 0; j < nx; ++j) 
					X[i][j] -= X[k][j] * m_factor[i][k];
							
		// Solve U*X = Y;
		for (int k = m_cols-1; k >= 0; --k) //int due to >= 0 condition
		{ 
			for (size_t j = 0; j < nx; ++j) 
			{
				X[k][j] /= m_factor[k][k];
			}
			for (size_t i = 0; i < size_t(k); ++i) 
				for (size_t j = 0; j < nx; ++j) 
					X[i][j] -= X[k][j] * m_factor[i][k];
		}
		return X;
	}

	stdVec_d LU::permute_copy(const stdVec_d &b) const
	{
		size_t piv_length = m_pivot.size();
		stdVec_d x(piv_length);
		for (size_t i = 0; i < piv_length; i++)
			x[i] = b[m_pivot[i]];
		return x;
	}

	mat LU::permute_copy(const mat &B) const
	{
		size_t piv_length = m_pivot.size();
		size_t nx = B[0].size();
		mat X(piv_length, nx);
		for (size_t i = 0; i < piv_length; ++i)
			for (size_t j = 0; j < nx; ++j)
				X[i][j] = B[m_pivot[i]][j];
		return X;
	}

	// -------------------------------------------------------
	//	Band Matrix Implementation
	// -------------------------------------------------------

	//constructor
	Band_matrix::Band_matrix(	size_t dim,
								size_t nu,
								size_t nl) :
								matrix(nu + nl + 2, dim),
								m_num_of_super(nu), //this gives the row index of the diagonal
								m_num_of_sub(nl),
								m_dimensions(dim),
								m_last_row_idx(nu + nl + 1),
								m_is_lu_decomposed(false)
	{}

	//Read access operator
	double Band_matrix::operator() (int i, int j) const {
		int row_idx = m_num_of_super - j + i;
		assert(row_idx >= 0 && row_idx < int(m_last_row_idx));
		return (*this)[size_t(row_idx)][size_t(j)];
	}

	// Write access operator
	double& Band_matrix::operator() (int i, int j) {
		int row_idx = m_num_of_super - j + i;
		assert(row_idx >= 0 && row_idx < int(m_last_row_idx));
		return (*this)[size_t(row_idx)][size_t(j)];
	}

	//	access to elements of the inverse of the diagonal (used in the decomposition) 
	//	saved as last row in band matrix
	//  Read value
	double Band_matrix::saved_diag(int j) const {
		assert(j < m_dimensions);
		return (*this)[m_last_row_idx][j];
	}
	//	Write value
	double & Band_matrix::saved_diag(int j) {
		assert(j < m_dimensions);
		return (*this)[m_last_row_idx][j];
	}

	// LR-Decomposition of a band matrix
	void Band_matrix::lu_decompose() {
		m_is_lu_decomposed = true;
		
		int  i_max, j_max;
		int  j_min;
		double x;
		// preconditioning
		// normalize the matrix such that a_ii = 1.0
		for (int i = 0; i < m_dimensions; ++i) 
		{
			assert(this->operator()(i, i) != 0.0); //TODO: find a better way to deal with diagonal elements == zero.
			saved_diag(i) = 1.0 / this->operator()(i, i);
			j_min = std::max(0, i - m_num_of_sub);
			j_max = std::min(m_dimensions - 1, i + m_num_of_super);

			for (int j = j_min; j <= j_max; ++j)
				this->operator()(i, j) *= saved_diag(i);
			
			this->operator()(i, i) = 1.0;          // avoids rounding errors
		}

		// Gauss LR-Decomposition
		for (int k = 0; k < m_dimensions; k++) 
		{
			i_max = std::min(m_dimensions - 1, k + m_num_of_sub);
			for (int i = k + 1; i <= i_max; i++) 
			{
				assert(this->operator()(k, k) != 0.0); //TODO: find a better way to deal with diagonal elements == zero.
				x = -this->operator()(i, k) / this->operator()(k, k);
				this->operator()(i, k) = -x;            // assembly part of L
				j_max = std::min(m_dimensions - 1, k + m_num_of_super);

				for (int j = k + 1; j <= j_max; ++j) // assembly part of R
					this->operator()(i, j) = this->operator()(i, j) + x*this->operator()(k, j);					
			}
		}
	}
	// solves Ly=b
	std::vector<double> Band_matrix::l_solve(const std::vector<double>& b) const 
	{
		assert(static_cast<size_t>(m_dimensions) == b.size());
		std::vector<double> x(m_dimensions);
		int j_start;
		double sum;
		for (int i = 0; i < m_dimensions; ++i) 
		{
			sum = 0;
			j_start = std::max(0, i - m_num_of_sub);

			for (int j = j_start; j < i; j++) 
				sum += this->operator()(i, j) * x[j];

			x[i] = (b[i] * saved_diag(i)) - sum;
		}
		return x;
	}
	// solves Rx=y
	std::vector<double> Band_matrix::r_solve(const std::vector<double>& b) const 
	{
		assert(static_cast<size_t>(m_dimensions) == b.size());
		std::vector<double> x(m_dimensions);
		int j_stop;
		double sum;
		for (int i = m_dimensions - 1; i >= 0; --i) //int due to >= 0 condition
		{ 
			sum = 0;
			j_stop = std::min(m_dimensions - 1, i + m_num_of_super);
			for (int j = i + 1; j <= j_stop; ++j) 
				sum += this->operator()(i, j) * x[j];
			x[i] = (b[i] - sum) / this->operator()(i, i);
		}
		return x;
	}

	std::vector<double> Band_matrix::lu_solve(const std::vector<double>& b) 
	{
		assert(static_cast<size_t>(m_dimensions) == b.size());
		std::vector<double>  x, y;
		if (m_is_lu_decomposed == false)
			lu_decompose();

		y = l_solve(b);
		x = r_solve(y);
		return x;
	}
}
