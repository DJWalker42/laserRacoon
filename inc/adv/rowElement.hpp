#ifndef ROWELEMENT_HPP
#define ROWELEMENT_HPP

#include <iostream>

namespace phys{

	template<class T> 
	class rowElement{
	public:
		rowElement(const T& val = 0, int col = -1);
		rowElement(const rowElement&e);
		~rowElement(){} // destructor

		const rowElement& operator=(const rowElement&e);

		const T& get_value() const
		{
			return value;
		} // read the value

		int get_column() const
		{
			return column;
		} // return the column

		const rowElement& operator+=(const T&t);

		const rowElement& operator+=(const rowElement<T>&e);

		const rowElement& operator-=(const T&t);

		const rowElement& operator-=(const rowElement<T>&e);

		const rowElement& operator*=(const T&t);

		const rowElement& operator/=(const T&t);
	private:
		T value;
		int column;
	};

	template<class T>
	bool operator<(const rowElement<T>&e, const rowElement<T>&f);

	template<class T>
	bool operator>(const rowElement<T>&e, const rowElement<T>&f);

	template<class T>
	bool operator==(const rowElement<T>&e, const rowElement<T>&f);

	template<class T>
	const rowElement<T> operator+(const rowElement<T>&e, const T&t);

	template<class T>
	const rowElement<T> operator+(const T&t, const rowElement<T>&e);

	template<class T>
	const rowElement<T> operator-(const rowElement<T>&e, const T&t);

	template<class T>
	const rowElement<T> operator*(const rowElement<T>&e, const T&t);

	template<class T>
	const rowElement<T> operator*(const T&t, const rowElement<T>&e);

	template<class T>
	const rowElement<T> operator/(const rowElement<T>&e, const T&t);

	template<class T>
	std::ostream& operator<<(std::ostream& os, const rowElement<T>& e);
}

#include <adv/rowElement.inl>

#endif