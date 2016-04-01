
namespace phys{

	template<class T>
	mesh<T>::mesh() : connectedList<T>(){}

	template<class T>
	mesh<T>::mesh(T& e) : connectedList<T>(e){} 

	template <class T>
	size_t mesh<T>::indexing()
	{
		for (mesh<T>* runner = this; runner; runner = (mesh<T>*) runner->next)
			runner->item.reset_indicies();
		size_t count = 0;
		for (mesh<T>* runner = this; runner; runner = (mesh<T>*) runner->next)
			runner->item.indexing(count);
		return count;
	}

}