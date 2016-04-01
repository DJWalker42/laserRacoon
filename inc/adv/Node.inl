namespace phys{

	template<class T>
	node<T>::node(	const T&loc,
					int idx,
					int sharing) :
					location(loc),
					index(idx),
					sharing_elements(sharing)
	{}

	template<class T>
	const node<T>& node<T>::operator=(const node<T>&n){
		if (this != &n){
			location = n.location;
			index = n.index;
			sharing_elements = n.sharing_elements;
		}
		return *this;
	} // assignment operator

	template<class T>
	std::ostream& operator<<(std::ostream& os, const node<T>& n)
	{
		os << n() << " index = " << n.get_index() << "; "
			<< n.get_sharing() << "\n";
		return os;
	}
}