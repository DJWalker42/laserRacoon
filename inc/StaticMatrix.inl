namespace phys{

	/************************* Constructors **************************************/
	template<class T, size_t M, size_t N>
	staticMat<T,M,N>::staticMat() : staticVec<staticVec<T, N>, M>(), m_trans(false) 
	{}

	template<class T, size_t M, size_t N>
	staticMat<T, M, N>::staticMat(	const staticVec<T, N>& v, 
									bool t) :
		staticVec<staticVec<T, N>, M>(v), m_trans(t) 
	{}

	template<class T, size_t M, size_t N>
	staticMat<T, M, N>::staticMat(	const staticVec<T, N>& v,
									const staticVec<T, N>& u,
									bool t ) :
		staticVec<staticVec<T, N>, M>(), m_trans(t)
	{
		(*this)[0] = v;
		(*this)[1] = u;
	}

	template<class T, size_t M, size_t N>
	staticMat<T, M, N>::staticMat(	const staticVec<T, N>& v,
									const staticVec<T, N>& u,
									const staticVec<T, N>& w,
									bool t) :
		staticVec<staticVec<T, N>, M>(), m_trans(t)
	{
		(*this)[0] = v;
		(*this)[1] = u;
		(*this)[2] = w;
	}

	template<class T, size_t M, size_t N>
	staticMat<T, M, N>::staticMat(	const T& val,
									bool t ) :
		staticVec<staticVec<T, N>, M>(), m_trans(t)
	{
		for (size_t i = 0; i < M; ++i)
			(*this)[i][i] = val;  						
	}

	template<class T, size_t M, size_t N>
	staticMat<T, M, N>::staticMat(	const staticVec<staticVec<T, N>, M>& vv,
									bool t ) :
		staticVec<staticVec<T, N>, M>(vv), m_trans(t)
	{}

	template<class T, size_t M, size_t N>
	staticMat<T, M, N>::staticMat(	const staticMat& m ) :
		staticVec<staticVec<T, N>, M>(m), m_trans(m.is_trans())
	{}

	/**************	Member functions and operators *********************/

	template<typename T, size_t M, size_t N>
	const staticMat<T, M, N>& staticMat<T, M, N>::operator=(const staticMat<T, M, N>& rhs)
	{
		if (this != &rhs)
		{
			*this = staticMat(rhs); //calls copy constructor is this wise?
		}
		return *this;
	}

	template<class T, size_t M, size_t N>
	const T& staticMat<T, M, N>::operator()(size_t i, size_t j) const
	{
		assert(i < M && j < N);
		return (*this)[i][j];
	}

	template<class T, size_t M, size_t N>
	T& staticMat<T, M, N>::operator()(size_t i, size_t j)
	{
		assert(i < M && j < N);
		return (*this)[i][j];
	}


	template<class T, size_t M, size_t N>
	const staticMat<T, M, N>& staticMat<T, M, N>::operator=(const T& val)
	{
		for (size_t i = 0; i < M; ++i)
			(*this)[i][i] = val;
		return *this;
	}

	template<class T, size_t M, size_t N>
	const staticMat<T, M, N>& staticMat<T, M, N>::operator*=(const staticMat<T, M, N>& m)
	{
		*this = matmul(*this, phys::transpose(m));
		return *this;
	}

	template<class T, size_t M, size_t N>
	const staticMat<T, M, N>& staticMat<T, M, N>::operator*=(const T& a)
	{
		for (size_t i = 0; i < M; ++i)
			(*this)[i] *= a;
		return *this;
	}

	template<class T, size_t M, size_t N>
	const staticMat<T, M, N>& staticMat<T, M, N>::operator/=(const T& a)
	{
		for (size_t i = 0; i < M; ++i)
			(*this)[i] /= a;
		return *this;
	}

	template<class T, size_t M, size_t N>
	const T staticMat<T, M, N>::l_one_norm(bool chk) const
	{
		T retval = 0;
		if (m_trans && !chk)
		{
			retval = this->l_inf_norm(m_trans);
		}
		else
		{
			for (size_t j = 0; j < N; ++j){
				T sum = 0;
				for (size_t i = 0; i < M; ++i){
					sum += fabs((*this)[i][j]);
				}
				if (sum > retval) retval = sum;
			}
		}
		return retval;
	}

	template<class T, size_t M, size_t N>
	const T staticMat<T, M, N>::l_inf_norm(bool chk) const
	{
		T retval = 0;
		if (m_trans && !chk)
		{
			retval = this->l_one_norm(m_trans);
		}
		else
		{
			for (size_t i = 0; i < M; ++i){
				T sum = 0;
				for (size_t j = 0; j < N; ++j){
					sum += fabs((*this)[i][j]);
				}
				if (sum > retval) retval = sum;
			}
		}
		return retval;
	}

	/******************* Non-member operators and functions *********************************/

	template<class T, size_t M, size_t N>
	const staticMat<T, M, N> operator*(const staticMat<T, M, N>& m, const T& a)
	{
		return staticMat<T, M, N>(m) *= a;
	}

	template<class T, size_t M, size_t N>
	const staticMat<T, M, N> operator*(const T& a, const staticMat<T, M, N>& m)
	{
		return staticMat<T, M, N>(m) *= a;
	}

	template<class T, size_t M, size_t N>
	const staticMat<T, M, N> operator/(const staticMat<T, M, N>& m, const T& a)
	{
		return staticMat<T, M, N>(m) /= a;
	}

	template<class T, size_t M, size_t N>
	const staticVec<T, N> operator*(const staticVec<T, M>& v, const staticMat<T, M, N>& m)
	{
		staticVec<T, N> result;
		for (size_t i = 0; i < M; ++i)
			result += v[i] * m[i];
		return result;
	}

	template<class T, size_t M, size_t N>
	const staticVec<T, M> operator*(const staticMat<T, M, N>& m,
		const staticVec<T, N>& v)
	{
		staticVec<T, M> result;
		for (size_t i = 0; i < M; ++i)
			result[i] = m[i] * v; //vector inner product
		return result;
	}

	template<class T, size_t M, size_t N, size_t K>
	const staticMat<T, M, K> operator*(const staticMat<T, M, N>&m1,
		const staticMat<T, N, K>&m2)
	{
		if (m1.is_trans())
		{
			if (!m2.is_trans())//m1.m_trans = true, m2.m_trans = false 
			{
				//this is the most computationally expensive command - try to avoid
				return matmul(transpose(m1), transpose(m2));
			}
			else//m1.m_trans = true, m2.m_trans = true
			{
				return matmul(transpose(m1), m2);
			}
		}
		else
		{
			if (!m2.is_trans())//m1.m_trans = false, m2.m_trans = false
			{
				return matmul(m1, transpose(m2));
			}
			else//m1.m_trans = false, m2.m_trans = true.
			{
				//this is the most efficient command -  try to use. 
				return matmul(m1, m2);
			}
		}
	}

	template<class T, size_t M, size_t N, size_t K>
	const staticMat<T, M, K> matmul(	const staticMat<T, M, N>& m1,
										const staticMat<T, K, N>& m2	) 
	{
		staticMat<T, M, K> result;
		for (size_t i = 0; i < M; ++i)
			result[i] = m2 * m1[i];
		return result;
	}

	template<class T, size_t M, size_t N, size_t K>
	void _matmul(staticMat<T, M, N>& m1, const staticMat<T, N, K>& m2)
	{
		for (size_t i = 0; i < M; ++i)
			m1[i] = m2 * m1[i]; //multiplication order important!
	}

	template<class T, size_t M, size_t N>
	const staticMat<T, N, M> transpose(const staticMat<T, M, N>&m)
	{
		//quick return if m empty or a single vector
		if (M < 2) return staticMat<T, M, N>(m);
		staticMat<T, N, M> result;
		//swap elements - here we have to scan through whole matrix
		for (size_t i = 0; i < M; ++i){
			for (size_t j = 0; j < N; ++j){
				result[j][i] = m[i][j];
			}
		}
		return result;
	}

	template<class T, size_t M, size_t N>
	void _transpose(staticMat<T, M, N>& m)
	{
		if (M < 2) return;
		assert(M == N); //Square matrices only; shape cannot change
		//swap elements
		T tmp;
		for (size_t i = 1; i < M; ++i){
			for (size_t j = 0; j < i; ++j){
				tmp = m[i][j];
				m[i][j] = m[j][i];
				m[j][i] = tmp;
			}
		}
	}

	template<class T, size_t M, size_t N>
	const T l_2_norm(const staticMat<T, M, N>& m)
	{
		return sqrt(m.l_one_norm() * m.l_inf_norm());
	}

	template<class T, size_t M, size_t N>
	std::ostream& operator<<(std::ostream& os, const staticMat<T, M, N>& m)
	{
		if (m.is_trans()){ //print out matrix as its transpose
			for (size_t j = 0; j < N; ++j){
				for (size_t i = 0; i < M; ++i){
					os << m[i][j];
				}
				os << "\n";
			}
		}
		else{
			for (size_t i = 0; i < M; ++i){
				os << m[i] << "\n";
			}
		}
		return os;
	}
}
