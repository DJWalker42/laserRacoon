	
namespace phys{

	/*********************************************************
	*	Constructors
	*********************************************************/
	template<class T, size_t N>
	staticVec<T, N>::staticVec( const T& a ) 
	{
		std::cout << "staticVec = " << a << std::endl;
		for (size_t i = 0; i < N; ++i)
			m_element[i] = a;
	}

	template<class T, size_t N>
	staticVec<T, N>::staticVec(	const T& a,
								const T& b )
	{
		std::cout << "staticVec double 2 constructor (" << a << "," << b << ")" << std::endl;
		m_element[0] = a;
		m_element[1] = b;
	}

	template<class T, size_t N>
	staticVec<T, N>::staticVec(	const T& a,
								const T& b,
								const T& c )
	{
		std::cout << "staticVec double 3 constructor (" 
				<< a << "," << b << "," << c << ")" << std::endl;
		m_element[0] = a;
		m_element[1] = b;
		m_element[2] = c;
	}


	template<class T, size_t N>
	staticVec<T, N>::staticVec(	const staticVec& v )
	{
		std::cout << "staticVec copy constructor" << std::endl;
		for (size_t i = 0; i < N; ++i)
			m_element[i] = v.m_element[i];
	}

	/***********************************************************
	*	Assignment operators
	***********************************************************/
	template<class T, size_t N>
	const staticVec<T, N>& staticVec<T, N>::operator=(const staticVec<T, N>& v)
	{
		if (this != &v)
			for (size_t i = 0; i < N; ++i)
				m_element[i] = v.m_element[i];
		return *this;
	}

	template<class T, size_t N>
	const staticVec<T, N>& staticVec<T, N>::operator=(const T& c)
	{
		for (size_t i = 0; i < N; ++i)
			m_element[i] = c;
		return *this;
	}

	/**********************************************************
	*	Access Operators
	**********************************************************/
	template<class T, size_t N>
	const T& staticVec<T, N>::operator[](size_t i) const
	{
		return m_element[i];
	} // read ith m_element of the vector

	template<class T, size_t N>
	T& staticVec<T, N>::operator[](size_t i)
	{
		return m_element[i];
	}// write ith m_element

	/*****************************************************************
	*	Member Operators
	******************************************************************/
	template<class T, size_t N>
	staticVec<T,N>& staticVec<T,N>::operator*=(const T& a)
	{
		for (size_t i = 0; i < N; ++i)
			m_element[i] *= a;
		return *this;
	}

	template<class T, size_t N>
	staticVec<T, N>& staticVec<T,N>::operator/=(const T& a)
	{
		std::cout << "staticVec /= T " << *this << std::endl;
		for (size_t i = 0; i < N; ++i)
			m_element[i] /= a;
		return *this;
	}

	template<class T, size_t N>
	staticVec<T, N>& staticVec<T,N>::operator+=(const staticVec& v)
	{
		std::cout << "staticVec += staticVec ";
		std::cout << *this << " + " << v << std::endl;
		for (size_t i = 0; i < N; ++ i){
			this->m_element[i] += v.m_element[i];
		}
		return *this;
	}

	template<class T, size_t N>
	staticVec<T,N>& staticVec<T,N>::operator-=(const staticVec& v)
	{
		for (size_t i = 0; i < N; ++ i){
			this->m_element[i] -= v.m_element[i];
		}
		return *this;
	}

	/*****************************************************************
	*	Non-member operators & functions
	******************************************************************/
	template<class T, size_t N>
	const staticVec<T, N> operator+(const staticVec<T, N>&u){
		return staticVec<T, N>(u);
	}//positive of a staticVec

	template<class T, size_t N>
	const staticVec<T, N> operator-(const staticVec<T, N>&u){
		return staticVec<T, N>(u) *= -1;
	} // negative of a staticVec

	template<class T, size_t N>
	const T operator*(const staticVec<T, N>&u, const staticVec<T, N>&v){
		T sum = 0;
		for (size_t i = 0; i < N; i++)
			sum += u[i] * +v[i];
		return sum;
	} // staticVec times staticVec (inner product)

	template<class T, size_t N>
	T squaredNorm(const staticVec<T, N>&u){
		return u*u;
	} // sum of squares of vector elements
}
