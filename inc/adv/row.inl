namespace phys{

	template<class T>
	row<T>::row() : connectedList<rowElement<T>>()
	{}

	template<class T>
	row<T>::row(const T&val, int col) : 
		connectedList<rowElement<T>>( rowElement<T>(val, col) )
	{}

	template<class T>
	row<T>::row(const connectedList<rowElement<T>>& c) :
		connectedList<rowElement<T>>(c)
	{}

	template<class T>
	const row<T>& row<T>::operator*=(const T&t)
	{
		item *= t;
		if (next) *(row*)next *= t;
		return *this;
	} // multiply by a T

	template<class T>
	const row<T>& row<T>::operator/=(const T&t)
	{
		item /= t;
		if (next) *(row*)next /= t;
		return *this;
	} // divide by a T

	template<class T>
	const T row<T>::operator*(const vector<T>& v) const
	{
		return next ? get_value() * v[get_column()] + *(row*)next * v
			: get_value() * v[get_column()];
	} // row times vector (inner product)

	template<class T>
	void row<T>::renumber_columns(const vector<int>& renumber)
	{
		item = rowElement<T>(get_value(), renumber[get_column()] - 1);
		if (next)
			(*(row<T>*)next).renumber_columns(renumber);
	} // renumber columns

	template<class T>
	const T row<T>::row_sum_coarse(const row<T>& r, const vector<int>& coarse) const
	{
		T contrib = coarse[get_column()] ? r[get_column()] : 0.;
		return next ? contrib + ((row<T>*) next)->row_sum_coarse(r, coarse) : contrib;
	}// row sum at coarse points

	template<class T>
	void row<T>::add_coarse(const row<T>& r, const vector<int>& coarse)
	{
		if (coarse[get_column()])
			item += r[get_column()];
		if (next) ((row<T>*) next)->add_coarse(r, coarse);
	}// row sum at coarse points

	template<class T>
	void row<T>::drop_items(const vector<int>& mask)
	{
		if (next)
		{
			if (!mask[(*next)().get_column()])
			{
				drop_next_item();
				drop_items(mask);
			}
			else
				(*(row<T>*)next).drop_items(mask);
			if (!mask[get_column()])drop_first_item();
		}
	} // "masking" the row by a vector of integers

	template<class T>
	void row<T>::drop_positive_items(	int i,
										const T& center,
										double thresh)
	{
		if (next)
		{
			if (((*next)().get_columns() != i) && ((*next)().get_value() / center >= -thresh))
			{
				drop_next_item();
				drop_positive_items(i, center, thresh);
			}
			else
				(*(row<T>*)next).drop_positive_items(i, center, thresh);

			if ((get_column() != i) && (get_value() / center >= -thresh))
				drop_first_item();
		}
	} // drop positive off-diagonal elements.

	template<class T>
	const row<T> operator*(const row<T>&r, const T&t)
	{
		return row<T>(r) *= t;
	} // row times T

	template<class T>
	const row<T> operator*(const T&t, const row<T>&r)
	{
		return row<T>(r) *= t;
	} // T times row

	template<class T>
	const row<T> operator/(const row<T>&r, const T&t)
	{
		return row<T>(r) /= t;
	} // row divided by a T

}