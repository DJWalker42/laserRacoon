#ifndef FINITE_ELEMENT_HPP
#define FINITE_ELEMENT_HPP

#include <adv/Node.hpp>

namespace phys{

	template<class T, size_t N>
	class FiniteElement{
	public:
		FiniteElement();

		FiniteElement(node<T>&, node<T>&, node<T>&);
		FiniteElement(FiniteElement<T,N>&);
		const FiniteElement<T,N>& operator=(FiniteElement<T,N>&);
		~FiniteElement(); //see definition below

		node<T>& operator()(size_t i){
			return *(vertex[i]);
		} // read/write ith vertex

		const node<T>& operator[](size_t i) const{
			return *(vertex[i]); 
		} // read only ith vertex 

		void reset_indicies(){
			for(size_t i = 0; i < N; ++i)
				vertex[i]->set_index(-1);
		}// reset indicies to -1

		void indexing( size_t& count )
		{
			for(size_t i = 0; i < N; ++i)
				if(vertex[i]->get_index() < 0)
					vertex[i]->set_index(count++);
		}
	private:
		node<T>* vertex[N];
	};

	// check whether a node n is in a finite element e ( n < e, reads "is n in e?")
	template<class T, size_t N>
	int operator<(const node<T>&n, const FiniteElement<T,N>&e);

	// print a finite element to screen, i.e. its nodes.
	template<class T, size_t N>
	std::ostream& operator<<(std::ostream& os, const FiniteElement<T,N>& e);

	typedef FiniteElement<point2, 3> triangle;		//	3 vertices in the cartesian plane (x,y)
	typedef FiniteElement<point3, 4> tetrahedron;	//	4 verticies in the cartesian space (x,y,z)
}

#include <adv/FiniteElement.inl>

#endif