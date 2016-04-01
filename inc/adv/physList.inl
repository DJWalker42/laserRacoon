namespace phys{

	/*****************************************************
	*	List class implementation
	*****************************************************/
	/* Constructors */

	template<class T>
	list<T>::list(size_t n) : 
		number(n), 
		item(n ? new T*[n] : 0)
	{}

	template<class T>
	list<T>::list(size_t n, const T&t) : 
		number(n), 
		item(n ? new T*[n] : 0)
	{
		for (size_t i = 0; i < number; ++i)
			item[i] = new T(t);
	}

	template<class T>
	list<T>::list(const list<T>&l) :
		number(l.number),
		item(l.number ? new T*[l.number] : 0)
	{
		for (size_t i = 0; i<l.number; i++)
			if (l.item[i]) item[i] = new T(*l.item[i]);
	}

	//assignment operator
	template<class T>
	const list<T>& list<T>::operator=(const list<T>& l)
	{
		if (this != &l){
			if (number > l.number)
				delete[](item + l.number);
			if (number < l.number){
				delete[] item;
				item = new T*[l.number];
			}

			for (size_t i = 0; i < l.number; i++)
				if (l.item[i]) item[i] = new T(*l.item[i]);
			number = l.number;
		}
		return *this;
	}

	//output stream window
	template<class T>
	std::ostream& operator<<(std::ostream& os, const list<T>& l)
	{
		for (size_t i = 0; i < l.size(); ++i)
			os << i << " = " << l[i] << "\n";
		return os;
	}

	/*****************************************************
	*	connectedList class implementation
	*****************************************************/
	/* Constructors */
	template<class T>
	connectedList<T>::connectedList() : 
		item(), 
		next(0) 
	{}

	template<class T>
	connectedList<T>::connectedList(T& t, connectedList* N) : 
		item(t), 
		next(N)
	{}

	template<class T>
	connectedList<T>::connectedList(const connectedList& l) :
		item(l()),
		next(l.next ? new connectedList(*l.next) : 0)
	{} 

	template<class T>
	connectedList<T>::~connectedList()
	{
		clear();
	}

	template<class T>
	const connectedList<T>&
		connectedList<T>::operator=(const connectedList<T>&L){
		if (this != &L){
			item = L();
			if (next){
				if (L.next)
					*next = *L.next;
				else{
					delete next;
					next = 0;
				}
			}
			else
				if (L.next) next = new connectedList(*L.next);
		}
		return *this;
	}

	template<class T>
	void connectedList<T>::drop_next_item(){
		if (next){
			if (next->next)
			{
				connectedList<T>* keep = next;
				next = next->next;
				keep->item.~T();
			}
			else
			{
				delete next;
				next = 0;
			}
		}
		else
			std::cout << "error: cannot drop nonexisting next item\n";
	}

	template<class T>
	void connectedList<T>::drop_first_item(){
		if (next){
			item = next->item;
			drop_next_item();
		}
		else
			std::cout << "error: cannot drop first item; no next.\n";
	}

	template<class T>
	void connectedList<T>::truncate_items(double threshold){
		if (next)
		{
			if (abs(next->item.get_value()) <= threshold)
			{
				drop_next_item();
				truncate_items(threshold);
			}
			else
				next->truncate_items(threshold);
		}

		if (next && (abs(item.get_value()) <= threshold))
			drop_first_item();
	}

	template<class T>
	const connectedList<T>& connectedList<T>::operator+=(connectedList<T>&L)
	{
		connectedList<T>* runner = this;
		connectedList<T>* Lrunner = &L;

		if (L.item < item){
			insert_first_item(L.item);
			Lrunner = L.next;
		}
		for (; runner->next; runner = runner->next){
			if (Lrunner && (Lrunner->item == runner->item)){
				runner->item += Lrunner->item;
				Lrunner = Lrunner->next;
			}
			for (; Lrunner && (Lrunner->item < runner->next->item);
				Lrunner = Lrunner->next){
				runner->insert_next_item(Lrunner->item);
				runner = runner->next;
			}
		}
		if (Lrunner && (Lrunner->item == runner->item)){
			runner->item += Lrunner->item;
			Lrunner = Lrunner->next;
		}
		if (Lrunner) runner->next = new connectedList<T>(*Lrunner);
		return *this;
	}

	template<class T>
	connectedList<T>& connectedList<T>::order(size_t length){
		if (length>1){
			connectedList<T>* runner = this;
			for (int i = 0; i < length / 2 - 1; i++)
				runner = runner->next;
			connectedList<T>* second = runner->next;
			runner->next = 0;
			order(length / 2);
			*this += second->order(length - length / 2);
		}
		return *this;
	}

	template<class T>
	std::ostream& operator<<(std::ostream& os, const connectedList<T>& l)
	{
		os << "item: \n" << l() << "\n";
		if (l.read_next()) os << *l.read_next();
		return os;
	}
}