#ifndef MULTIGRID_HPP
#define MULTIGRID_HPP

#include <adv/SparseMat.h>

namespace phys{

	const double grid_ratio = 0.95; 
	const int useILU = 0;

	template<class T> class Multigrid;

	template<class T>
	std::ostream& operator<<(std::ostream& os, const Multigrid<T>& mg);

	template<class T>
	class Multigrid
	{
	public:
		Multigrid();
		Multigrid(const Multigrid& mg);
		Multigrid(const sparseMat<T>& m );

		~Multigrid();

		const Multigrid<T>& operator=(const Multigrid<T>&);

		const vector<T>& Vcycle (const vector<T>&, vector<T>&) const;
		friend std::ostream& operator<< <T>(std::ostream& os, const Multigrid<T>&);
		
	private:
		sparseMat<T> A;
		sparseMat<T> U;
		sparseMat<T> L;
		sparseMat<T> P;
		sparseMat<T> R;
		Multigrid* next;
	};

	/*	
		Nu1 is the number of relaxations to apply before the coarse-grid corrections. 
		Nu2 is the number of relaxations to apply after the coarse-grid corrections.
		cycle_index is the number of coarse-grid corrections to apply.
		Nu_coarse is the number of relaxations to apply at the coarsest grid.
	*/
	const int Nu1 = 1, Nu2 = 1, cycle_index = 1, Nu_coarse = 1;
}


#include <adv/Multigrid.inl>



#endif