namespace phys{

	template<class T>
	rowElement<T>::rowElement(	const T& val, 
								int col ) : 
								value(val), 
								column(col)
	{}//constructor

	template<class T>
	rowElement<T>::rowElement(	const rowElement&e): 
								value(e.value), 
								column(e.column)
	{} // copy constructor

	template<class T>
	const rowElement<T>& rowElement<T>::operator=(const rowElement&e)
	{
		if (this != &e)
		{
			value = e.value;
			column = e.column;
		}
		return *this;
	} // assignment operator

	template<class T>
	const rowElement<T>& rowElement<T>::operator+=(const T&t)
	{
		value += t;
		return *this;
	} // adding a T

	template<class T>
	const rowElement<T>& rowElement<T>::operator+=(const rowElement<T>&e)
	{
		value += e.value;
		return *this;
	} // adding a rowElement

	template<class T>
	const rowElement<T>& rowElement<T>::operator-=(const T&t)
	{
		value -= t;
		return *this;
	} // subtracting a T

	template<class T>
	const rowElement<T>&
		rowElement<T>::operator-=(const rowElement<T>&e)
	{
		value -= e.value;
		return *this;
	} // subtracting a rowElement

	template<class T>
	const rowElement<T>&
		rowElement<T>::operator*=(const T&t)
	{
		value *= t;
		return *this;
	} // multiplying by a T

	template<class T>
	const rowElement<T>& rowElement<T>::operator/=(const T&t)
	{
		value /= t;
		return *this;
	} // dividing by a T


	/* Non-member operators */
	template<class T>
	bool operator<(const rowElement<T>&e, const rowElement<T>&f)
	{
		return e.get_column() < f.get_column();
	} // smaller column index

	template<class T>
	bool operator>(const rowElement<T>&e, const rowElement<T>&f)
	{
		return e.get_column() > f.get_column();
	} // greater column index

	template<class T>
	bool operator==(const rowElement<T>&e, const rowElement<T>&f)
	{
		return e.get_column() == f.get_column();
	} // same column

	template<class T>
	const rowElement<T> operator+(const rowElement<T>&e, const T&t)
	{
		return rowElement<T>(e) += t;
	} // rowElement plus a T

	template<class T>
	const rowElement<T> operator+(const T&t, const rowElement<T>&e)
	{
		return rowElement<T>(e) += t;
	} // T plus rowElement

	template<class T>
	const rowElement<T> operator-(const rowElement<T>&e, const T&t)
	{
		return rowElement<T>(e) -= t;
	} // rowElement minus T

	template<class T>
	const rowElement<T> operator*(const rowElement<T>&e, const T&t)
	{
		return rowElement<T>(e) *= t;
	} // rowElement times a T

	template<class T>
	const rowElement<T> operator*(const T&t, const rowElement<T>&e)
	{
		return rowElement<T>(e) *= t;
	} // T times rowElement

	template<class T>
	const rowElement<T> operator/(const rowElement<T>&e, const T&t)
	{
		return rowElement<T>(e) /= e;
	} // rowElement divided by a T

	template<class T>
	std::ostream& operator<<(std::ostream& os, const rowElement<T>& e)
	{
		os << e.get_value() << " : column index = " << e.get_column();
		return os;
	}



}