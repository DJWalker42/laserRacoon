#ifndef PHYSLIST_HPP
#define PHYSLIST_HPP

#include <iostream>
#include <algorithm>

namespace phys{

	/**********************************************	
	/		list class 
	**********************************************/
	template<class T> class list{
	public:
		list(size_t n=0);

		list(size_t n, const T&t);

		list(const list<T>&);

		const list<T>& operator=(const list<T>&);

		~list(){
			for(size_t i = 0; i < number; ++i)
				delete item[i];
			delete [] item;
		} // destructor

		size_t size() const{
			return number;
		} // list size

		T& operator()(size_t i){
			if(item[i]) return *(item[i]);
			else throw std::runtime_error( "List index out-of-bounds" );
		} // read/write ith item

		const T& operator[](size_t i)const{
			if(item[i]) return *(item[i]);
			else throw std::runtime_error( "List index out-of-bounds" );
		} // read only ith item
	protected:
		size_t number;
		T** item;
	};

	//output stream operator for lists
	template<class T>
	std::ostream& operator<<(std::ostream& os, const list<T>& l);

	/**********************************************	
	/		connectedList class 
	**********************************************/
	template<class T> 
	class connectedList{
	public:
		connectedList();

		connectedList(T& t, connectedList* N=0);

		connectedList(	const connectedList& l );

		~connectedList();

		const connectedList& operator=(const connectedList&);

		const T& operator()() const{ return item;}

		const connectedList* read_next() const{return next;}

		connectedList& last(){
			return next ? next->last() : *this;
		} // last item

		const size_t length() const{
			return next ? next->length() + 1 : 1;
		} // number of items

		void append(T& t){
			last().next = new connectedList(t);
		} // append item

		void insert_next_item(T& t){
			next = new connectedList(t,next);
		} // insert item in second place

		void insert_first_item(T& t){
			next = new connectedList(item,next);
			item = t;
		} // insert item at the beginning

		void clear()
		{
			while(next)
				drop_next_item();
		}// clears list up to first item

		void drop_next_item();
		void drop_first_item();
		void truncate_items(double);
		const connectedList& operator+=(connectedList&);
		connectedList& order(size_t);
	protected:
		T item;
		connectedList* next;
	};

	template<class T>
	std::ostream& operator<<(std::ostream& os, const connectedList<T>& l);
}

#include <adv/physList.inl>

#endif