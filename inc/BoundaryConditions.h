#ifndef DJW_BOUNDARYCONDITIONS_H
#define DJW_BOUNDARYCONDITIONS_H

namespace phys{

	//forward declaration of the BvpODE class
	class BvpODE;

	enum boundaryType {DIRICHLET, NEUMANN};//, CAUCHY}; //note to implement Cauchy BCs we'd need _TWO_ doubles for each boundary (value _AND_ derviative).

	/*	BoundaryConditions: a class to store the type and value of boundary conditions use when solving boundary value problems in one-dimension. As we are
		in one dimension the boundary values are necessarly singluar values. 
		
		For boundary value problems in 2 dimensions the boundary conditions are, in general, functions along the specific dimension. 
		We would have to store these as function pointers as data members in the class, one function pointer for each boundary (top, bottom, left, right).
		Each boudary also needs to specify its type. 
	*/
	class BoundaryConditions{
		friend class BvpODE;
	public:	
		BoundaryConditions(boundaryType lhsType, double lhsValue, boundaryType rhsType, double rhsValue) :
			m_LHStype(lhsType),
			m_RHStype(rhsType),
			m_LHSvalue(lhsValue),
			m_RHSvalue(rhsValue)
		{}

		/* Use to change the boundary condition on the left-hand-side of the domain */
		void setLhsCondition(boundaryType bcType, double bcValue)
		{
			m_LHStype = bcType;
			m_LHSvalue = bcValue;
		}

		/* Use to change the boundary condition on the right-hand-side of the domain */
		void setRhsCondition(boundaryType bcType, double bcValue)
		{
			m_RHStype = bcType;
			m_RHSvalue = bcValue;
		}

	private:
		boundaryType m_LHStype;
		boundaryType m_RHStype;

		double m_LHSvalue;
		double m_RHSvalue;
	};

}

#endif
