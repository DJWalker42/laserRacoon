namespace phys{

	template<class T>
	Multigrid<T>::Multigrid() : 
		A(),
		U(),
		L(),
		P(),
		R(),
		next(0)
	{}

	template<class T>
	Multigrid<T>::Multigrid(const Multigrid& mg): 
							A(mg.A),
							U(mg.U),
							L(mg.L),
							P(mg.P),
							R(mg.R),
							next(mg.next ? new Multigrid(*mg.next) : 0)
	{}

	template<class T>
	Multigrid<T>::Multigrid(const sparseMat<T>& m) : 
		A(m),
		U(useILU ? A : 0),
		L(useILU ? U.factorize(0.05) : 0),
		P(A),
		R(P.create_transfer()),
		next(R.row_number() <= grid_ratio*A.row_number() ?
				new Multigrid(R*A*P) : 0)
	{}

	template<class T>
	Multigrid<T>::~Multigrid()
	{
		delete next;
		next = 0;
	}

	template<class T>
	const Multigrid<T>& Multigrid<T>::operator=(const Multigrid<T>& mg)
	{
		if (this != &mg)
		{
			A = mg.A;
			U = mg.U;
			L = mg.L;
			P = mg.P;
			R = mg.R;
			if (next) //target has a coarser grid
			{
				if (mg.next) //operand has a coarser grid
					*next = *mg.next; //recursive call
				else // remove coarser grid
				{
					delete next;
					next = 0;
				}
			}
			else // add new coarse grid if operand has coarser grid.
				if (mg.next) next = new Multigrid(*mg.next);
		}
		return *this;
	}// asignment operator

	template<class T>
	std::ostream& operator<<(std::ostream& os, const Multigrid<T>& mg)
	{
		os << "A =\n" << mg.A << "\n";
		os << "P =\n" << mg.P << "\n";
		os << "R =\n" << mg.R << "\n";
		if (mg.next) os << "\n" << *mg.next; //recursive call
	}

	template<class T>
	const vector<T>& Multigrid<T>::Vcycle (	const vector<T>& f,
											vector<T>& x) const
	{
		if (next) // apply recursion to the next Multigrid
		{
			//do the fisrt Nu1 relaxations
			for (int i = 0; i < Nu1; ++i)
			{
				if (useILU)
					ILU(A, L, U, f, x);
				else
					Gauss_Seidel(A, L, U, f, x);
			}

			vector<T> residual = f - A * x;
			vector<T> correction(R.row_number());
			//apply coarse-grid corrections
			for (int i = 0; i < cycle_index; ++i)
				next->Vcycle(R * residual, correction);
			x += P * correction;
			//apply Nu2 relaxations
			for (int i = 0; i < Nu2; ++i)
			{
				if (useILU)
					ILU(A, L, U, f, x);
				else
					reverse_Gauss_Seidel(A, L, U, f, x);
			}
		}
		else //at the coarsest level
		{
			for (int i = 0; i < Nu_coarse; ++i)
			{
				if (useILU)
					ILU(A, L, U, f, x);
				else
					symmetric_Gauss_Seidel(A, L, U, f, x);
			}
		}
	}

}