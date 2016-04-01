namespace phys{

	/* Constructors */

	template<class T>
	sparseMat<T>::sparseMat(size_t n) :
		list<row<T>>(n)
	{
		for (size_t i = 0; i < n; ++i)
			item[i] = 0;
	}

	template<class T>
	sparseMat<T>::sparseMat(size_t n, const T& a) :
		list<row<T>>(n)
	{
		for (size_t i = 0; i < n; ++i)
			item[i] = new row<T>(a, i);
	}

	template<class T>
	sparseMat<T>::sparseMat(mesh<triangle>&m)
	{
		item = new row<T>*[number = m.indexing()];
		for (size_t i = 0; i < number; ++i)
			item[i] = 0;

		point2 gradient[3];
		gradient[0] = point2(-1, -1);
		gradient[1] = point2(1, 0);
		gradient[2] = point2(0, 1);

		for (const mesh<triangle>* runner = &m; runner; runner = (const mesh<triangle>*)runner->read_next())
		{
			/*	(*runner)()		gets the current triangle.
			[i]				gets the ith node of that triangle.
			()				gets the location (x,y) of that node --or--
			.get_index()	the index of that node
			*/

			/*	Code interprets the arguments of the staticMat constructor as rows of the matrix,
			thus we set the final argument to true to indicate that we wish to treat
			S as its transpose. */
			staticMat2 S((*runner)()[1]() - (*runner)()[0](),
				(*runner)()[2]() - (*runner)()[0](), true);
			staticMat2 S_inv = inverse(S);
			//due to the way matmul works matmul(A,A) returns with the result A * A'
			staticMat2 weight = abs(det(S) / 2.) * matmul(S_inv, S_inv);
			for (size_t i = 0; i < 3; ++i)
				for (size_t j = 0; j < 3; ++j)
				{

					int I = (*runner)()[i].get_index();
					int J = (*runner)()[j].get_index();
					if (item[I])
					{
						row<T> r(gradient[j] * weight*gradient[i], J);
						*item[I] += r;
					}
					else
						item[I] = new row<T>(gradient[j] * weight*gradient[i], J);
				}
		}
	}

	/****************** Member functions ************************/

	template<class T>
	int sparseMat<T>::column_number() const
	{
		int maxColumn = -1;
		for (int i = 0; i < row_number(); ++i)
			if (item[i]) maxColumn = std::max(maxColumn, item[i]->last()().get_column());
		return maxColumn + 1;
	} // number of columns

	template<class T>
	const sparseMat<T>& sparseMat<T>::operator+=(const sparseMat<T>& M)
	{
		for (int i = 0; i < row_number(); ++i)
			*item[i] += *M.item[i];
		return *this;
	} // add a sparse matrix to current sparse matrix

	template<class T>
	const sparseMat<T>& sparseMat<T>::operator-=(const sparseMat<T>&M)
	{
		for (int i = 0; i < row_number(); i++)
		{
			row<T> minus = -1. * *M.item[i];
			*item[i] += minus;
		}
		return *this;
	} // subtract a sparse matrix from current sparse matrix

	template<class T>
	const sparseMat<T>& sparseMat<T>::operator*=(const T& t)
	{
		for (int i = 0; i < row_number(); i++)
			*item[i] *= t;
		return *this;
	} // multiply by T

	template<class T>
	const sparseMat<T>& sparseMat<T>::truncate(double t)
	{
		int rows = this->row_number();
		for (int i = 0; i < rows; ++i){
			this->item[i]->truncate_items(t * abs((*this)(i, i)));
		}
		return *this;
	}//remove off-diagonal elements in rows .LE. to t * current row diagonal element.

	template<class T>
	const sparseMat<T> sparseMat<T>::factorize(double thresh)
	{
		sparseMat<T> L(row_number());
		L.item[0] = new row<T>(0., 0);

		for (int i = 0; i < row_number(); ++i)
		{
			while (item[i] && (item[i]->get_column() < i))
			{
				T factor = item[i]->get_value() / item[item[i]->get_column()]->get_value();
				if (abs(factor) >= thresh)
				{
					*item[i] += (-factor) * *item[item[i]->get_column()];
					if (L.item[i])
						L.item[i]->append(factor, item[i]->get_column());
					else
						L.item[i] = new row<T>(factor, item[i]->get_column());
				}

				item[i]->drop_first_item();
			}

			if (!L.item[i]) L.item[i] = new row<T>(0., 0);
			this->truncate(thresh);
		}
		return L;
	}// incomplete LU factorisation. Note: this will change current matrix to incomplete U, returning incomplete L

	template<class T>
	const vector<T> sparseMat<T>::forward_elimination(const vector<T>& f) const
	{
		vector<T> result(f.size());
		result[0] = f[0];
		for (size_t i = 1; i < f.size(); ++i)
			result[i] = item[i] ? f[i] - *item[i] * result : f[i];
		return result;
	}// forward elimination in L

	template<class T>
	const vector<T> sparseMat<T>::back_substitution(const vector<T>& f) const
	{
		vector<T> result(f.size(), 0.);
		for (int i = f.size() - 1; i >= 0; --i)
		{
			result[i] = item[i]->read_next() ? f[i] - *(row<T>*)item[i]->read_next() * result : f[i];
			result[i] /= (*item[i])().get_value();
		}
		return result;
	}// backward substitution in U

	template<class T>
	const vector<int> sparseMat<T>::coarsen(double thresh) const
	{
		/*	This constructs the coarse grid, c, for the multigrid method.
		The coarser grid is a subset of the finer grid, and it follows
		the rule that the ith component of the finer grid only has a
		non-zero index in the coarser grid if and only if i is in c.
		*/
		//initialise c to contain all the finer grid points
		vector<int> coarse(row_number(), 1);
		//drop points that don't meet threshold conditions
		for (int i = 0; i < row_number(); ++i)
		{
			if (coarse[i])
				for (const row<T>* runner = item[i]; runner; runner = (const row<T>*)runner->read_next())
					if ((runner->get_column() != i) && (abs(runner->get_value()) >= thresh * abs((*this)(i, i))))
						coarse[runner->get_column()] = 0;
		}
		//picks up on opposite parity values of off-diagonal to diagonal elements, and adds back in to c if necessary
		for (int i = 0; i < row_number(); ++i)
		{
			if (!coarse[i])
			{
				int drop = 1;
				for (const row<T>* runner = item[i]; runner; runner = (const row<T>*)runner->read_next())
					if ((coarse[runner->get_column()]) && (runner->get_value() / (*this)(i, i) <= -thresh))
						drop = 0;
				if (drop) coarse[i] = 1;
			}
		}
		// for all i in c order the indicies in sequence.
		int count = 0;
		for (int i = 0; i < row_number(); ++i)
			if (course[i]) coarse[i] = ++count;
		return coarse;
	}// define the coarse grid


	template<class T>
	const sparseMat<T> sparseMat<T>::create_transfer()
	{
		const vector<int> coarse = coarsen();
		for (int i = 0; i < row_number(); ++i)
		{
			if (coarse[i])
			{
				delete item[i];
				item[i] = new row<T>(1., coarse[i] - 1);
			}
			else
			{
				item[i]->drop_postive_items(i, (*item[i])(i), 0.05);
				item[i]->drop_items(coarse);
				item[i]->renumber_columns(coarse);
				*item[i] /= item[i]->row_sum();
			}
		}
		return transpose(*this);
	}// create transfer operators

	/****************** Non-member functions ******************************/

	template<class T>
	const sparseMat<T> operator+(const sparseMat<T>& M1, const sparseMat<T>& M2)
	{
		return sparseMat<T>(M1) += M2;
	} 

	template<class T>
	const sparseMat<T> operator-(const sparseMat<T>& M1, const sparseMat<T>& M2)
	{
		return sparseMat<T>(M1) -= M2;
	} 

	template<class T>
	const sparseMat<T> operator*(const T& t, const sparseMat<T>& M)
	{
		return sparseMat<T>(M) *= t;
	}

	template<class T>
	const sparseMat<T> operator*(const sparseMat<T>&M, const T&t)
	{
		return sparseMat<T>(M) *= t;
	} 

	template<class T>
	const vector<T> operator*(const sparseMat<T>&M, const vector<T>&v)
	{
		vector<T> result(M.row_number(), 0.);
		for (int i = 0; i < M.row_number(); ++i)
			result(i) = M[i] * v;
		return result;
	} 

	template<class T>
	const sparseMat<T> operator*(const sparseMat<T>& m1, const sparseMat<T>& m2)
	{
		sparseMat<T> result(m1.row_number());
		for (int i = 0; i < m1.row_number(); ++i)
		{
			result.item[i] = new row<T>(m1.item[i]->get_value() *
				*m2.item[m1.item[i]->get_column()]);
			for (const row<T>* runner = (const row<T>*) m1.item[i]->read_next();
				runner; runner = (const row<T>*) runner->read_next())
			{
				row<T> r = runner->get_value() * *m2.item[runner->get_column()];
				*result.item[i] += r;
			}
		}
		return result;
	}

	template<class T>
	void Gauss_Seidel(const sparseMat<T>& A,
		const vector<T>& f,
		vector<T>& x)
	{
		for (size_t i = 0; i < f.size(); ++i)
			x[i] += (f[i] - *A.item[i] * x) / A(i, i);
	}

	template<class T>
	void reverse_Gauss_Seidel(const sparseMat<T>& A,
		const vector<T>& f,
		vector<T>& x)
	{
		for (int i = int(f.size() - 1); i >= 0; --i) //int due to >= 0 condition
			x[i] += (f[i] - *A.item[i] * x) / A(i, i);
	}

	template<class T>
	void symmetric_Gauss_Seidel(const sparseMat<T>& A,
		const vector<T>& f,
		vector<T>& x)
	{
		vector<T> diagA(f.size());
		for (size_t i = 0; i < f.size(); ++i){
			diagA[i] = A(i, i);
			x[i] += (f[i] - *A.item[i] * x) / diagA[i];
		}
		for (int i = int(f.size() - 2); i >= 0; --i) //f.size()-2 correct
			x[i] += (f[i] - *A.item[i] * x) / diagA[i];
	}

	template<class T>
	const vector<T> operator/(const vector<T>& v,
		const sparseMat<T>& A)
	{
		vector<T> result(v);
		for (size_t i = 0; i < v.size(); ++i)
			result[i] /= A(i, i); // divide by zero??
		return result;
	}

	template<class T>
	void Jacobi(const sparseMat<T>& A,
		const vector<T>& f,
		vector<T>& x)
	{
		for (size_t i = 0; i < f.size(); ++i)
			x += (f - A*x) / A;
	}

	template<class T>
	const sparseMat<T> transpose(const sparseMat<T>&M)
	{
		sparseMat<T> Mt(M.column_number());
		for (int i = 0; i < M.row_number(); i++)
			for (const row<T>* runner = M.item[i]; runner; runner = (const row<T>*)runner->read_next())
			{
				if (Mt.item[runner->get_column()])
					Mt.item[runner->get_column()]->append(runner->get_value(), i);
				else
					Mt.item[runner->get_column()] = new row<T>(runner->get_value(), i);
			}
		return Mt;
	} 

	template<class T>
	const sparseMat<T> diagonal(const sparseMat<T>&M)
	{
		size_t rows = M.row_number();
		sparseMat<T> diag(rows);
		for (size_t i = 0; i < rows; ++i)
			diag.item[i] = new row<T>(M(i, i), i);
		return diag;
	}

	template<class T>
	void ILU(const sparseMat<T>& A,
		const sparseMat<T>& L,
		const sparseMat<T>& U,
		const vector<T>& f,
		vector<T>& x)

	{
		x += U.back_substitution(L.forward_elimination(f - A * x));
	} 

	template<class T>
	std::ostream& operator<<(std::ostream& os, const sparseMat<T>& m)
	{

		for (int i = 0; i < m.row_number(); ++i){
			for (int j = 0; j < m.column_number(); ++j){
				os << std::setprecision(3);
				os << std::setfill(' ');
				os << std::setw(7);
				os << std::fixed;
				os << m(i, j); //if element doesn't exist this returns a zero
			}
			os << "\n";
		}
		return os;
	}

}