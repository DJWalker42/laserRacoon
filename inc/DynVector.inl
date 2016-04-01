namespace phys{
	/************************************************************************
	*	Implementation of overloaded operators and others
	************************************************************************/
	template<typename T>
	std::vector<T>& operator+= (std::vector<T>& u, const std::vector<T>& v)
	{
		assert(u.size() = v.size());
		size_t n = u.size();
		for(size_t i = 0; i < n; ++i)
			u[i] += v[i]; 
		return u;
	}

	template<typename T>
	std::vector<T>& operator-= (std::vector<T>& u, const std::vector<T>& v)
	{
		assert(u.size() = v.size());
		size_t n = u.size();
		for (size_t i = 0; i < n; ++i)
			u[i] -= v[i];
		return u;
	}

	template<typename T>
	std::vector<T>& operator+=(std::vector<T>& v, const T& a)
	{
		size_t n = v.size();
		for(size_t i = 0; i < n; ++i)
			v[i] += a;
		return v;
	}

	template<typename T>
	std::vector<T>& operator-=(std::vector<T>& v, const T& a)
	{
		size_t n = v.size();
		for (size_t i = 0; i < n; ++i)
			v[i] -= a;
		return v;
	}

	template<typename T>
	std::vector<T>& operator*= (std::vector<T>& u, const T& a)
	{
		size_t n = u.size();
		for(size_t i = 0; i < n; ++i)
			u[i] *= a; 
		return u;
	}

	template<typename T>
	std::vector<T>& operator/= (std::vector<T>& u, const T& a)
	{
		size_t n = u.size();
		for (size_t i = 0; i < n; ++i)
			u[i] /= a;
		return u;
	}

	template<typename T>
	const std::vector<T> operator+(const std::vector<T>&u, const std::vector<T>&v)
	{
		return std::vector<T>(u) += v;
	}

	template<typename T>
	const std::vector<T> operator-(const std::vector<T>&u, const std::vector<T>&v)
	{
		return std::vector<T>(u) -= v;
	}

	template<typename T>
	const std::vector<T> operator+(const std::vector<T>&u)
	{
		return std::vector<T>(u);
	}

	template<typename T>
	const std::vector<T> operator-(const std::vector<T>&u)
	{
		return std::vector<T>(u) *= -1;
	}

	template<typename T>
	const std::vector<T> operator*(const std::vector<T>&u, const T& a)
	{
		return std::vector<T>(u) *= a;
	}

	template<typename T>
	const std::vector<T> operator*(const T& a, const std::vector<T>&u)
	{
		return std::vector<T>(u) *= a;
	}

	template<typename T>
	const std::vector<T> operator/(const std::vector<T>& u, const T& a)
	{
		return std::vector<T>(u) /= a;
	}

	template<typename T>
	const T operator*(const std::vector<T>&u, const std::vector<T>&v)
	{
		assert(u.size() == v.size());
		size_t dim = u.size();
		T sum = 0;
		for (size_t i = 0; i < dim; i++)
			sum += u[i] * +v[i];
		return sum;
	}

	template<typename T>
	const T squaredNorm(const std::vector<T>&u)
	{
		return u*u;
	}

	template<typename T>
	const T norm(const std::vector<T>&v, int p)
	{
		T sum = 0;
		for(int i = 0; i < v.size(); ++i)
		{
			sum += pow(fabs(v[i]), p);
		}
		return pow(sum, T(1)/T(p)); 
	}

	template<typename T>
	std::ostream & operator<<(std::ostream& os, const std::vector<T>& v)
	{
		if (v.empty()) 
			os << "NULL\n";
		else
			for (size_t i = 0; i < v.size(); ++i)
				os << v[i] << " ";
		return os;
	}

	template<typename T>
	const std::vector<T> read_vector(const std::string& line, const T& t)
	{
		std::vector<T> retval;
		import_vector(line, retval);
		return retval;
	}

	template<typename T>
	const std::vector<T> sub_vector( const std::vector<T>& v, size_t start, size_t end)
	{
		assert(end < v.size() &&  start <= end ); 
		typename std::vector<T>::const_iterator first = v.begin() + start;
		typename std::vector<T>::const_iterator last = v.begin() + end + 1;

		return std::vector<T>(first, last); 
	}
}
