#ifndef NODE_HPP
#define NODE_HPP

#include <staticVector.hpp>

namespace phys{

	/* T will typically be a point2 or point3 depending on the application */
	template<class T> class node{
	public:
		node(const T&loc = 0., int idx = -1, int sharing = 0);

		const node& operator=(const node&);

		~node(){}//destructor

		const T& operator()()const{
			return location;
		}//get "coordinate" location

		int get_index() const{
			return index;
		}

		void set_index(int i){
			index = i;
		}

		int get_sharing() const{
			return sharing_elements;
		}

		void more_sharing(){
			++sharing_elements;
		}

		bool less_sharing(){
			return !(sharing_elements ? --sharing_elements : 0);
		}

		bool singleton() const{
			return !sharing_elements;
		}// indicate a dangling node - ohh matron!!
	private:
		T location;
		int index;
		int sharing_elements;
	};

	template<class T>
	std::ostream& operator<<(std::ostream& os, const node<T>& n);
}

#include <adv/Node.inl>

#endif