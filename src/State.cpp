#include "State.h"

#include <stdexcept>

namespace phys{
	namespace ode{

		/* Default constructor*/
		state::state() : x(), y(), dy() {}

		/*constructor for 1st order system with one dimension*/
		state::state(double independent_var, double dependent_var)
			: x(independent_var)
		{
			y.push_back(dependent_var);
			num_dims = 1;
		}

		/*constructor for 2nd order system with one dimension*/
		state::state(double independent_var, double dependent_var, double first_deriv)
			: x(independent_var)
		{
			y.push_back(dependent_var);
			dy.push_back(first_deriv); //used for checking state
			y.push_back(first_deriv); //data actually used here.
			num_dims = 1; 
		}

		/*constructor for 1st or 2nd order system with multiple dimensions*/
		state::state(double independent_var, stdVec_d dependent_vars, stdVec_d first_deriv)
			: x(independent_var), y(dependent_vars), dy(first_deriv)
		{
			if (dy.size() != y.size()) __error_message_dimensions();
			num_dims = y.size();
			y.insert(y.end(), dy.begin(), dy.end());
		}

		void state::__error_message_dimensions()
		{
			std::string errmsg = "Error: mismatch between the number of variables\n";
			errmsg += "in the dependent vector and the first derivative vector.\n";
			errmsg += "They should have the same size.";
			throw(std::runtime_error(errmsg));
		}


	}
}
