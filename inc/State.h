#ifndef STATE_HPP
#define STATE_HPP

#include <string>

#include "DynVector.h"

namespace phys{
	namespace ode{

		/**	Structure for the state of a system under the application of a first or second ordered ODE.
			Multiple dimensions can be handled. The member vector y will contain both the dependent variable(s)
			and, when appropriate, the first derivative varaible(s) after construction. The member vector dy is
			not modified by the ODE_solver class; the first derivative values are handled within the vector y.
		*/
		struct state{
			//representation
			double x;				//!<independent variable
			stdVec_d y;				//!<dependent variable(s)
			stdVec_d dy;			//!<first derivative(s) of the dependent variable(s)
			unsigned int num_dims;	//!<number of dimensions to for which to solve. Auto set by input size of y.

			/* Default constructor*/
			state();

			/*constructor for 1st order system with one dimension*/
			state(double independent_var, double dependent_var);

			/*constructor for 2nd order system with one dimension*/
			state(double independent_var, double dependent_var, double first_deriv);

			/*constructor for 1st or 2nd order system with multiple dimensions*/
			state(double independent_var, stdVec_d dependent_vars,
				stdVec_d first_deriv = stdVec_d());

		private:
			void __error_message_dimensions();
		};
	}
}

#endif
