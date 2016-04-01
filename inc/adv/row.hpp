#ifndef ROW_HPP
#define ROW_HPP

#include <iostream>
#include <adv/physList.hpp>
#include <physVector.hpp>
#include <adv/rowElement.hpp>

namespace phys{

	template<class T>
	class row : public connectedList<rowElement<T> >{
	public:
		row();

		row(const T&, int);

		row(const connectedList<rowElement<T>>&);

		const rowElement<T>& operator()() const
		{
			return item;
		} // read only first item

		const T& get_value() const
		{
			return item.get_value();
		} // read only first-item value

		int get_column() const
		{
			return item.get_column();
		} // read only first-item column

		void insert_next_item(const T&val, int col)
		{
			rowElement<T> e(val,col);
			connectedList<rowElement<T> >::insert_next_item(e);
		} // insert a rowElement as second item

		void insert_first_item(const T&val, int col)
		{
			rowElement<T> e(val,col);
			connectedList<rowElement<T> >::insert_first_item(e);
		} // insert a rowElement at the beginning

		void append(const T&val, int col){
			rowElement<T> e(val,col);
			connectedList<rowElement<T> >::append(e);
		} // append a rowElement at the end of row

		const T row_sum() const
		{
			return next ? get_value() + (*(row<T>*)next).row_sum() : get_value();
		} // row-sum

		const T operator[](int i) const
		{
			return (get_column() == i) ? get_value() : next && (get_column() < i) ? (*(row*)next)[i] : 0.;
		} // read only the value at column i

		const row& operator*=(const T&t);
		const row& operator/=(const T&t);

		const T operator*(const vector<T>& v) const;

		void renumber_columns(const vector<int>& renumber);

		const T row_sum_coarse( const row<T>& r, const vector<int>& coarse ) const;

		void add_coarse( const row<T>& r, const vector<int>& coarse);

		void drop_items(const vector<int>&);
		void drop_positive_items(int, const T&, double);
	};

	template<class T>
	const row<T> operator*(const row<T>&r, const T&t);

	template<class T>
	const row<T> operator*(const T&t, const row<T>&r);

	template<class T>
	const row<T> operator/(const row<T>&r, const T&t);
}

#include <adv/row.inl>

#endif