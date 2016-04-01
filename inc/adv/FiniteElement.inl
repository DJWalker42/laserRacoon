namespace phys{

	template<class T, size_t N>
	FiniteElement<T,N>::FiniteElement()
	{
		for (size_t i = 0; i < N; ++i)
			vertex[i] = new node<T>;
	} 

	template<class T, size_t N>
	FiniteElement<T, N>::FiniteElement(	node<T>&a,
										node<T>&b,
										node<T>&c )
	{
		vertex[0] = a.singleton() ? new node<T>(a) : &a;
		vertex[1] = b.singleton() ? new node<T>(b) : &b;
		vertex[2] = c.singleton() ? new node<T>(c) : &c;
		for (size_t i = 0; i < N; ++i){
			vertex[i]->more_sharing();
		}
	}

	template<class T, size_t N>
	FiniteElement<T, N>::FiniteElement(FiniteElement<T, N>& e)
	{
		for (size_t i = 0; i < N; ++i)
		{
			vertex[i] = e.vertex[i];
			vertex[i]->more_sharing();
		}
	}

	template<class T, size_t N>
	const FiniteElement<T, N>& FiniteElement<T, N>::operator=(FiniteElement<T, N>& e)
	{
		if (this != &e){
			//delete vertices not sharing elements first.
			for (size_t i = 0; i < N; ++i)
				if (vertex[i]->less_sharing())
					delete vertex[i];
			for (size_t i = 0; i < N; ++i){
				vertex[i] = e.vertex[i];
				vertex[i]->more_sharing();
			}
		}
		return *this;
	}

	template<class T, size_t N>
	FiniteElement<T, N>::~FiniteElement(){
		for (size_t i = 0; i < N; ++i)
			if (vertex[i]->less_sharing())
				delete vertex[i];
	}// destructor

	template<class T, size_t N>
	int operator<(const node<T>&n, const FiniteElement<T, N>&e){
		for (size_t i = 0; i < N; ++i)
			if (&n == &(e[i]))
				return int(i + 1);
		return 0;
	}// @todo: check if we change return type to bool

	template<class T, size_t N>
	std::ostream& operator<<(std::ostream& os, const FiniteElement<T, N>& e)
	{
		for (size_t i = 0; i < N; ++i)
			os << e[i];
		return os;
	}
}