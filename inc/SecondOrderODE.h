#ifndef DJW_SECONDORDERODE_H
#define DJW_SECONDORDERODE_H

namespace phys{

	//forward declaration of the BvpODE class
	class BvpODE;

	/*	
		Class to describe a second order ode in one dimension -	
		- Can this class also describe a first-order ode? If so, how?
		- How might we extend this to cover higher dimensions?	
		- How will you deal with the crossed derivative i.e. Uxy? 
		Ignoring it could be a valid answer, depends on the context; do crossed derivatives occur often in physical systems?	
	*/
	class SecondOrderODE{
		friend class BvpODE; 
	public:	
		SecondOrderODE(	double Uxx, 
						double Ux, 
						double U, 
						double (*pDiff) (double), 
						double xMin, 
						double xMax) :
						m_Uxx(Uxx),
						m_Ux(Ux),
						m_U(U),
						m_pDiff(pDiff),
						m_xMin(xMin),
						m_xMax(xMax)
		{}

	private:
		double m_Uxx;					//!< Coefficient of Uxx
		double m_Ux;					//!< Coefficient of Ux
		double m_U;						//!< Coefficient of U

		double (*m_pDiff) (double x);	//!< function pointer to the differential (rhs) function (note the single argument for 1D) 
		double m_xMin;					//!< The lhs boundary of the solution domain
		double m_xMax;					//!< The rhs boundary of the solution domain
	};

}

#endif