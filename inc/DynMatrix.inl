namespace phys{

	/*************** Constructors ***********************/
	template<class T>
	matrix<T>::matrix( size_t m ) : data_array( m, std::vector<T>(m, 0) ), trans(false)
	{}

	template<class T>
	matrix<T>::matrix(size_t m, size_t n, bool is_T) :
		data_array(m, std::vector<T>(n)), trans(is_T) 
	{}

	template<class T>
	matrix<T>::matrix(size_t m, size_t n, const T& val, bool is_T) :
		data_array(m, std::vector<T>(n, val)), trans(is_T) 
	{}

	template<class T>
	matrix<T>::matrix(size_t m, const std::vector<T>& row, bool is_T) :
		data_array(m, row), trans(is_T)
	{}

	template<class T>
	matrix<T>::matrix(const std::vector<std::vector<T>>& vv, bool is_T) :
		data_array(vv), trans(is_T)
	{}

	template<class T>
	matrix<T>::matrix(const std::vector<T>& v, size_t rows, size_t cols, bool is_T) :
		data_array(rows, std::vector<T>(cols)),
		trans(is_T)
	{
		assert(v.size() == rows*cols);
		if (trans)
		{
			for (size_t i = 0; i < rows; ++i)
				for (size_t j = 0; j < cols; ++j)
					data_array[i][j] = v[j*rows + i];
			trans = false; // we have read it in as transposed 
		}
		else
		{
			for (size_t i = 0; i < rows; ++i)
				for (size_t j = 0; j < cols; ++j)
					data_array[i][j] = v[i*cols + j];
		}

	}

	template<class T>
	matrix<T>::matrix(const matrix& mat) :
		data_array(mat.data_array), trans(mat.trans)
	{}


	template<typename T>
	matrix<T>& matrix<T>::operator=( const matrix& m )
	{
		if(this != &m)
		{
			data_array = m.data_array;
			trans = m.trans;
		}
		return *this;
	}

	/************** Member functions and operators **********************/

	template<class T>
	const T& matrix<T>::operator()(size_t i, size_t j) const
	{
		return data_array[i][j]; 
	}

	template<typename T>
	T& matrix<T>::operator()(size_t i, size_t j)
	{
		return data_array[i][j]; 
	}

	template<typename T>
	const std::vector<T>& matrix<T>::operator[](size_t i) const
	{
		return data_array[i];
	}

	template<typename T>
	std::vector<T>& matrix<T>::operator[](size_t i)
	{
		return data_array[i];
	}

	template<class T>
	const matrix<T> matrix<T>::operator()(size_t imin, size_t imax, size_t jmin, size_t jmax) const
	{
		assert(imax >= imin && jmax >= jmin);
		size_t irange = (imax - imin) + 1;
		size_t jrange = (jmax - jmin) + 1;
		matrix<T> result(irange, jrange);
		for (size_t i = 0; i < irange; ++i){
			for (size_t j = 0; j < jrange; ++j){
				result[i][j] = data_array[imin + i][jmin + j];
			}
		}
		return result;
	}

	template<class T>
	const matrix<T>& matrix<T>::operator=(const T& val)
	{
		size_t imax = data_array.size();
		if (val == 1) //set up the identity matrix (ones in main diagonal)
			for (size_t i = 0; i < imax; ++i)
				data_array[i][i] = val;
		else//broadcast value to all elements
			for (size_t i = 0; i < imax; ++i)
				data_array[i] = val;
		return *this;
	}

	template<typename T>
	const matrix<T>& matrix<T>::operator+=(const T& a)
	{
		for(size_t i = 0; i < data_array.size(); ++i)
			data_array[i] += a;
		return *this;
	}

	template<typename T>
	const matrix<T>& matrix<T>::operator-=(const T& a)
	{
		for (size_t i = 0; i < data_array.size(); ++i)
			data_array[i] -= a;
		return *this;
	}

	template<typename T>
	const matrix<T>& matrix<T>::operator+=(const matrix<T>& m)
	{

		assert(data_array.size() == m.data_array.size());
		for (size_t i = 0; i < data_array.size(); ++i)
			data_array[i] += m.data_array[i]; 
		return *this;
	}

	template<typename T>
	const matrix<T>& matrix<T>::operator-=(const matrix<T>& m)
	{
		assert(data_array.size() == m.data_array.size());
		for (size_t i = 0; i < data_array.size(); ++i)
			data_array[i] -= m.data_array[i];
		return *this;
	}

	//matrix member operator functions
	template<class T>
	const matrix<T>& matrix<T>::operator*=(const T& a)
	{
		size_t imax = data_array.size();
		for (size_t i = 0; i < imax; ++i)
			data_array[i] *= a;
		return *this;
	}//multiplication by scalar

	template<class T>
	const matrix<T>& matrix<T>::operator/=(const T& a)
	{
		size_t imax = data_array.size();
		for (size_t i = 0; i < imax; ++i)
			data_array[i] /= a;
		return *this;
	}//division by scalar

	template<class T>
	const matrix<T>& matrix<T>::swap_rows(size_t i, size_t j, bool chk)
	{
		if (trans && !chk)
		{
			this->swap_cols(i, j, trans);
		}
		else
		{
			if (i != j)
			{
				std::vector<T> keep = data_array[i];
				data_array[i] = data_array[j];
				data_array[j] = keep;
			}
		}
		return *this;
	}

	template<class T>
	const matrix<T>& matrix<T>::swap_cols(size_t i, size_t j, bool chk)
	{
		if (trans && !chk)
		{
			this->swap_rows(i, j, trans);
		}
		else
		{
			if (i != j)
			{
				size_t rows = this->size();
				for (size_t k = 0; k < rows; ++k)
				{
					T keep = data_array[k][i];
					data_array[k][i] = data_array[k][j];
					data_array[k][j] = keep;
				}
			}
		}
		return *this;
	}

	template<typename T>
	const matrix<T>& matrix<T>::resize(size_t resize_size, const std::vector<T>& val, dimension h_w)
	{
		switch(h_w)
		{
		case HEIGHT:
			if(!data_array.empty() && !val.empty())
				assert( val.size() == data_array[0].size() );
			//if empty we just append the row.
			data_array.resize(resize_size, val); 
			break;
		case WIDTH:
			if(data_array.empty())
			{
				data_array.resize(val.size());//append empty row vectors
			}
			else
				assert( val.size() == data_array.size() );

			size_t rows = data_array.size();
			for (size_t i = 0; i < rows; ++i)
				data_array[i].resize(resize_size, val[i]);

			break;
		}
		return *this;
	}

	template<class T>
	const matrix<T>& matrix<T>::pad_matrix(size_t rows, size_t cols, const T& val)
	{
		size_t trow = (trans) ? cols : rows;
		size_t tcol = (trans) ? rows : cols;

		//resize number of rows (this either adds empty row vectors or removes current ones)
		data_array.resize(trow);

		//then resize number of columns
		for (size_t i = 0; i < trow; ++i)
			data_array[i].resize(tcol, val); //this either removes columns or adds new ones with val.

		return *this;
	}

	template<typename T>
	void matrix<T>::append_row( const std::vector<T>& v )
	{
		//check dimensions
		if(!data_array.empty())
			assert( data_array[0].size() == v.size() );
		data_array.push_back(v); 
	}

	template<typename T>
	void matrix<T>::append_col(const std::vector<T>& v)
	{
		//check dimensions
		if (!data_array.empty())
			assert(data_array.size() == v.size());
		for(size_t i = 0; i < data_array.size(); ++i)
			data_array[i].push_back( v[i] );
	}

	template<typename T>
	void matrix<T>::clear()
	{
		data_array.clear();
	}

	template<typename T>
	bool matrix<T>::empty()
	{
		return data_array.empty();
	}
	
	/*********************************************************
		Non-member Functions 	
	*********************************************************/
	template<typename T>
	const matrix<T> operator+(const matrix<T>&m1, const matrix<T>&m2)
	{
		return matrix<T>(m1) += m2;
	}

	template<typename T>
	const matrix<T> operator-(const matrix<T>&m1, const matrix<T>&m2)
	{
		return matrix<T>(m1) -= m2;
	}

	template<typename T>
	const matrix<T> operator+(const T&a, const matrix<T>&m)
	{
		return matrix<T>(m) += a;
	}

	template<typename T>
	const matrix<T> operator+(const matrix<T>&m, const T&a)
	{
		return matrix<T>(m) += a;
	}

	template<typename T>
	const matrix<T> operator-(const T&a, const matrix<T>&m)
	{
		return -matrix<T>(m) += a;
	}

	template<typename T>
	const matrix<T> operator-(const matrix<T>&m, const T&a)
	{
		return matrix<T>(m) -= a; 
	}

	template<typename T>
	const matrix<T> operator+(const matrix<T>&m)
	{
		return matrix<T>(m); 
	}

	template<typename T>
	const matrix<T> operator-(const matrix<T>&m)
	{
		return matrix<T>(m) *= -1;
	}

	template<class T>
	const matrix<T> operator*(const T&a, const matrix<T>&m)
	{
		return matrix<T>(m) *= a;
	}

	template<class T>
	const matrix<T> operator*(const matrix<T>&m, const T&a)
	{
		return matrix<T>(m) *= a;
	}

	template<class T>
	const matrix<T> operator/(const matrix<T>&m, const T&a)
	{
		return matrix<T>(m) /= a;
	}

	template<class T>
	const std::vector<T> operator*(const std::vector<T>&v, const matrix<T>&m)
	{
		assert(v.size() == m.size());
		size_t imax = m.size();
		std::vector<T> result(m[0].size());
		for (size_t i = 0; i < imax; ++i)
			result += v[i] * m[i]; //scalar times vector returns a vector
		return result;
	}

	template<class T>
	const std::vector<T> operator*(const matrix<T>&m, const std::vector<T>&v)
	{
		assert(v.size() == m[0].size());
		size_t imax = m.size();
		std::vector<T> result(imax);
		for (size_t i = 0; i < imax; ++i)
			result[i] = m[i] * v; //vector inner product returns a T
		return result;
	}

	template<class T>
	const matrix<T> operator*(const matrix<T>&m1, const matrix<T>&m2)
	{
		/*
		assert(m1[0].size() == m2.size());
		size_t imax = m1.size();
		matrix<T> result(imax, m2[0].size());
		for(size_t i = 0; i < imax; ++i)
			result[i] = m1[i] * m2;
		return result;
		*/

		if (m1.is_trans())
		{
			if (!m2.is_trans())//m1.trans = true, m2.trans = false 
			{
				//this is the most computationally expensive call - try to avoid
				return matmul(transpose(m1), transpose(m2));
			}
			else//m1.trans = true, m2.trans = true
			{
				return matmul(transpose(m1), m2);
			}
		}
		else
		{
			if (!m2.is_trans())//m1.trans = false, m2.trans = false
			{
				return matmul(m1, transpose(m2));
			}
			else//m1.trans = false, m2.trans = true.
			{
				//this is the most efficient call -  try to use. 
				return matmul(m1, m2);
			}
		}
	}

	template<class T>
	const matrix<T> matmul(const matrix<T>& m1, const matrix<T>& m2)
	{
		assert(m1[0].size() == m2[0].size());//note that m2 is the TRANSPOSE of the actual rhs operand matrix
		size_t imax = m1.size();
		matrix<T> result(imax, m2.size());
		for (size_t i = 0; i < imax; ++i)
			result[i] = m2 * m1[i];
		return result;
	}

	template<class T>
	void _matmul(matrix<T>& m1, const matrix<T>& m2)
	{
		assert(m1[0].size() == m2[0].size());
		size_t imax = m1.size();
		for (size_t i = 0; i < imax; ++i)
			m1[i] = m2 * m1[i];
	}

	template<class T>
	const matrix<T> transpose(const matrix<T>&m)
	{
		//quick return if m empty or a single vector
		if (m.size() < 2) return matrix<T>(m);

		size_t M = m.size(), N = m[0].size();
		matrix<T> result(N, M);
		//assign elements to new positions
		for (size_t i = 0; i < M; ++i){
			for (size_t j = 0; j < N; ++j){
				result[j][i] = m[i][j];
			}
		}
		return result;
	}

	template<class T>
	void _transpose(matrix<T>& m)
	{
		//quick return if m empty or a single vector
		if (m.size() < 2) return;

		size_t M = m.size(), N = m[0].size();
		size_t min_dim = phys::maths::min(M, N);

		//first do the square section of the matrix
		//We only need to scan through strictly lower triangle 
		//hence start at second row (i = 1) and avoid i == j (diagonal)
		for (size_t i = 1; i < min_dim; ++i){
			for (size_t j = 0; j < i; ++j){
				T tmp = m[i][j];
				m[i][j] = m[j][i];
				m[j][i] = tmp;
			}
		}
		//then deal with "extra" rectangular section
		if (M < N)
			// height < width: landscape -> portrait
			// columns in extra width become rows of extra height.
		{
			for (size_t j = min_dim; j < N; ++j){
				m.append_row(std::vector<T>()); //append an empty row to m; will be the m[j]th row
				for (size_t i = 0; i < min_dim; ++i){
					m[j].push_back(m[i][j]);  //fill the new jth row with values from the jth column
				}
			}
			//remove extra width (columns) from each row in original m
			for (size_t i = 0; i < min_dim; ++i){
				m[i].resize(min_dim);
			}
		}
		else if (M > N)
			// height > width: portrait -> landscape
			// rows in extra height become columns of extra width
		{
			for (size_t j = 0; j < N; ++j){
				for (size_t i = min_dim; i < M; ++i){
					m[j].push_back(m[i][j]); // append values from the jth column to the jth row.
				}
			}
			//remove extra rows from original m
			m.resize(min_dim);
		} //else square -- all done.
		//all that remains to be done is to toggle the trans flag such that if m was previously being treated as its transpose
		//then it should be toggled to false as m is now physically its transpose, else do nothing (trans flag already false).
		if (m.is_trans()) m.transpose();
	}

	template<class T>
	const matrix<T> identity(size_t dim)
	{
		matrix<T> I(dim); //square matrix filled with zeros
		for(size_t i = 0; i < dim; ++i)
			I(i,i) = 1; //assign one to the diagonal elements
		return I; 
	}

	template<class T>
	std::ostream& operator<<(std::ostream& os, const matrix<T>& m)
	{
		if (m.empty()) { os << "NULL\n"; return os; }

		if (m.is_trans()){ //print out matrix as its transpose
			for (size_t j = 0; j < m[0].size(); ++j){
				for (size_t i = 0; i < m.size(); ++i){
					os << m[i][j] << " ";
				}
				os << "\n";
			}
		}
		else{
			for (size_t i = 0; i < m.size(); ++i){
				os << m[i] << "\n";
			}
		}
		return os;
	}

	template<class T>
	const matrix<T> swap_rows(const matrix<T>&m, size_t i, size_t j)
	{
		return matrix<T>(m).swap_rows(i, j);
	}

	template<class T>
	const matrix<T> swap_cols(const matrix<T>&m, size_t i, size_t j)
	{
		return matrix<T>(m).swap_cols(i, j);
	}

	template<class T>
	const matrix<T> resize(const matrix<T>&m, size_t resize_size, const std::vector<T>& val, dimension which)
	{
		return matrix<T>(m).resize_rows(resize_size, val, which);
	}

	template<class T>
	const matrix<T> pad_matrix(const matrix<T>&m, size_t rows, size_t cols, const T& val)
	{
		return matrix<T>(m).pad_matrix(rows, cols, val);
	}

	template<class T>
	matrix<T> read_matrix(const std::string& filename, const T& t, bool transpose)
	{
		size_t rows = 0, cols = 0;
		std::vector<T> temp;
		import_array(filename, temp, rows, cols);
		if (transpose) swap(rows, cols);
		return matrix<T>(temp, rows, cols, transpose);
	}
}
