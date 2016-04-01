namespace phys{

	template<class T>
	Polynomial<T>::Polynomial(size_t n) : 
		list<T>(n)
	{
		for (size_t i = 0; i < n; i++)
			item[i] = 0;
	}

	template<class T>
	Polynomial<T>::Polynomial(size_t n, const T&a) :
		list<T>(n, a)
	{}

	template<class T>
	const Polynomial<T>&operator+=(Polynomial<T>& p, const Polynomial<T>&q)
	{
		if (p.size() >= q.size())
			for (size_t i = 0; i < q.size(); i++)
				p(i) += q[i];
		else{
			Polynomial<T> keepQ = q;
			p = keepQ += p;
		}
		return p;
	} 

	template<class T>
	const Polynomial<T>operator+(const Polynomial<T>& p,
		const Polynomial<T>&q)
	{
		Polynomial<T> keep = p;
		return keep += q;
	} 

	template<class T>
	const Polynomial<T>&operator*=(Polynomial<T>& p, const T&a)
	{
		for (size_t i = 0; i < p.size(); i++)
			p(i) *= a;
		return p;
	} 

	template<class T>
	const Polynomial<T>operator*(const T&a, const Polynomial<T>&p)
	{
		Polynomial<T> keep = p;
		return keep *= a;
	} 

	template<class T>
	Polynomial<T>operator*(const Polynomial<T>&p, const Polynomial<T>&q)
	{
		Polynomial<T> result(p.degree() + q.degree() + 1, 0);
		for (size_t i = 0; i < result.size(); i++)
			for (size_t j = std::max(0, int(i) - int(q.degree()));
				j <= std::min(i, p.degree()); j++)
		{
			if (j == std::max(0, int(i) - int(q.degree())))
				result(i) = p[j] * q[i - j];
			else
				result(i) += p[j] * q[i - j];
		}
		return result;
	} 

	template<class T>
	Polynomial<T>&operator*=(Polynomial<T>&p, const Polynomial<T>&q)
	{
		return p = p * q;
	} 

	template<class T>
	const T calculatePolynomial(const Polynomial<T>&p, const T&x)
	{
		T powerOfX = 1;
		T sum = 0;
		for (size_t i = 0; i < p.size(); i++){
			sum += p[i] * powerOfX;
			powerOfX *= x;
		}
		return sum;
	} 

	template<class T>
	const T HornerPolynomial(const Polynomial<T>&p, const T&x)
	{
		T result = p[p.degree()];
		for (size_t i = p.degree(); i > 0; --i)
		{
			result *= x;
			result += p[i - 1];
		}
		return result;
	} 

	template<class T>
	const list<T>deriveRinverse(const T&r, size_t n)
	{
		list<T> Rinverse(n + 1, 0);
		Rinverse(0) = 1 / r;
		for (size_t i = 0; i<n; i++)
			Rinverse(i + 1) = -double(i + 1) / r * Rinverse[i];
		return Rinverse;
	}

	template<class T>
	const T Taylor(const list<T>&f, const T&h)
	{
		T pow_of_h = 1;
		T sum = 0;
		for (size_t i = 0; i < f.size() - 1; ++i){
			sum += f[i] * pow_of_h;
			pow_of_h *= h / (i + 1);
		}
		return sum;
	} 

	template<class T>
	const T HornerTaylor(const list<T>&f, const T&h)
	{
		T result = f[f.size() - 2];
		for (size_t i = f.size() - 2; i > 0; --i){
			result *= h / i;
			result += f[i - 1];
		}
		return result;
	} 

	template<class T>
	const T deriveProduct(const list<T>&f, const list<T>&g, size_t n)
	{
		T sum = 0;
		if (n < Pascal::n_max){ //use the pre-calculated values if we have them
			for (size_t i = 0; i <= n; i++)
				sum += Pascal::triangle[n - i].pos[i] * f[i] * g[n - i];
		}
		else{ //compute the binomial coefficients
			std::vector<int> pascal_triangle = Pascal::compute_triangle(n);
			for (size_t i = 0; i < n; i++)
				sum += pascal_triangle[i] * f[i] * g[n - i - 1];
		}
		return sum;
	} 

	template<class T>
	const T integral(const Polynomial<T>&p)
	{
		T sum = 0;
		for (size_t i = 0; i < p.size(); i++)
			sum += (1. / (i + 1)) * p[i];
		return sum;
	} 

	template<class T>
	const T integral(const Polynomial<Polynomial<T> >&p)
	{
		Polynomial<T> sum(p.size() + 1, 0);
		Polynomial<T> one(1, 1);
		Polynomial<T> x(2, 0);
		x(1) = 1;
		Polynomial<T> oneMinusX(2, 1);
		oneMinusX(1) = -1;
		list<Polynomial<T> > xPowers(p.size(), one);
		list<Polynomial<T> > oneMinusXpowers(p.size() + 1, one);
		for (int i = 1; i<p.size(); i++)
			xPowers(i) = x * xPowers[i - 1];
		for (int i = 1; i <= p.size(); i++)
			oneMinusXpowers(i) = oneMinusX * oneMinusXpowers[i - 1];
		for (int k = p.degree(); k >= 0; k--)
			for (int j = 0; j <= k; j++)
				sum += (p[k][j] / (j + 1))
				* oneMinusXpowers[j + 1] * xPowers[k - j];
		return integral(sum);
	} 
}