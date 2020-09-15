#include <exception>
#include <stdexcept>
#include <cmath>

#include "RootSearch.h"

namespace phys{
	namespace roots{

		/**************************************************************************
		*	RootSearch abstract Base Class
		**************************************************************************/
		//Constructors
		RootSearch::RootSearch() : m_function(), m_derivative(), m_lft_brc(), m_rht_brc()
		{ initialise_consts(); }

		RootSearch::RootSearch(double(*func)(double)) : m_function(func)
		{ initialise_consts(); }

		RootSearch::RootSearch(double(*func)(double), double left, double right) :
			m_function(func), m_lft_brc(left), m_rht_brc(right)
		{ initialise_consts(); }

		RootSearch::RootSearch(double(*func)(double), double(*deriv)(double)) :
			m_function(func), m_derivative(deriv) 
		{ initialise_consts(); }

		RootSearch::RootSearch(double(*func)(double), double(*deriv)(double), double initial_guess) :
			m_function(func), m_derivative(deriv), m_lft_brc(initial_guess), m_rht_brc(initial_guess)
		{ initialise_consts(); }

		RootSearch::RootSearch(double(*func)(double), double(*deriv)(double), double left, double right) :
			m_function(func), m_derivative(deriv), m_lft_brc(left), m_rht_brc(right)
		{ initialise_consts(); }

		void RootSearch::initialise_consts()
		{
			m_tolerance = 1.e-8;
			m_iteration_count = 0;
			m_iter_max = 100;
			m_error = 1.0;
			m_is_bisec = true;	
			m_is_newton = false;
		}

		//Utility functions
		bool RootSearch::find_brackets(double start, double step, double end)
		{
			if (start >= end) return false;
			check_function_set(__FUNCTION__);
			m_lft_brc = start;
			bool retval = false;
			do{
				m_rht_brc = m_lft_brc + step;
				if (m_function(m_lft_brc) * m_function(m_rht_brc) < 0){
					retval = true;
				} else {
					m_lft_brc = m_rht_brc;
				}
			} while (step*(m_lft_brc - end) < 0 && false == retval);
			return retval;
		}

		void RootSearch::set_brackets(double new_lft_brc, double new_rht_brc)
		{
			check_function_set(__FUNCTION__);
			m_lft_brc = new_lft_brc;
			m_rht_brc = new_rht_brc;
			if (m_is_bisec)
				check_root_is_bracketed(__FUNCTION__, m_function(m_lft_brc), m_function(m_rht_brc));
			//setting new brackets implies looking for a new root so reset error to 1
			m_error = 1.0;
		}

		void RootSearch::check_all(std::string call_func)
		{
			check_function_set(call_func);
			if (m_is_newton)
				check_derivative_set(call_func);
			check_brackets_set(call_func);
			if (m_is_bisec)
				check_root_is_bracketed(call_func, m_function(m_lft_brc), m_function(m_rht_brc));
		}

		void RootSearch::check_function_set(const std::string& calling_function)
		{
			if (m_function == 0)
			{
				std::string errmsg = "Function not set.\n";
				errmsg += "Use set_function(<func_name>) before trying\n";
				errmsg += "to use " + calling_function;
				throw (std::runtime_error(errmsg));
			}
		}

		void RootSearch::check_derivative_set(const std::string& calling_function)
		{
			if (m_derivative == 0)
			{
				std::string errmsg = "Derivative not set.\n";
				errmsg += "Use set_derivative(<deriv_name>) before trying\n";
				errmsg += "to use " + calling_function;
				throw (std::runtime_error(errmsg));
			}
		}

		void RootSearch::check_brackets_set(const std::string& call_func)
		{
			//FIXME: what if zero is a valid starting bracket value?
			if (m_lft_brc == double() || m_rht_brc == double())
			{
				std::string errmsg = call_func + ": Brackets have not been set.\n";
				errmsg += "Please either explicitly set_brackets(l,r) or\n";
				errmsg += " use find_bracket(start, step, end)";
				throw (std::runtime_error(errmsg));
			}
		}

		void RootSearch::check_root_is_bracketed(	const std::string& call_func,
													double func_lft,
													double func_rht	)
		{
			if (func_lft * func_rht > 0)
			{
				std::string errmsg = call_func + ": Root not bracketed";
				throw (std::runtime_error(errmsg));
			}
		}

		void RootSearch::output_data(const std::string& calling_func, double root)
		{
			std::cout << "Root search: " << calling_func << "\n";
			std::cout << "Root computed at: " << root << "\n";
			std::cout << "Error achieved: " << m_error << "\n";
			std::cout << "Iteration count: " << m_iteration_count << "\n\n";
		}


		double Bisection::find_root(bool print_output)
		{
			std::string func_name = "BISECTION";
			check_all(func_name);
			double func_lft = m_function(m_lft_brc);
			double root = (m_lft_brc + m_rht_brc)/2;
			double func_mid = m_function(root);
			size_t iter_cnt = 0;
			while(m_error > m_tolerance && iter_cnt++ < m_iter_max)
			{				
				if(func_lft * func_mid < 0)	//root contained in left subinterval
				{
					m_rht_brc = root;
				}else						//root contained in right subinterval
				{
					m_lft_brc = root;
					func_lft = func_mid;
				}			
				root = (m_lft_brc + m_rht_brc)/2;
				m_error = fabs((m_rht_brc - m_lft_brc)/2/root);
				func_mid = m_function(root); 
			}
			//store iter_count for output
			m_iteration_count = iter_cnt;
			if(print_output) output_data(func_name, root);
			return root;
		}

		double Secant::find_root(bool print_output)
		{
			std::string func_name = "SECANT";
			check_all(func_name);
			double root = m_rht_brc;
			double x1 = m_lft_brc;
			double f0 = m_function(x1), f1;
			double delta;
			size_t iter_cnt = 0; 
			while(m_error > m_tolerance && iter_cnt++ < m_iter_max)
			{
				f1 = m_function(root);
				delta = -f1 * (root - x1)/(f1 - f0);
				x1 = root;
				f0 = f1;
				root += delta;
				m_error = fabs(delta/root);
			}
			//store iter_count for output
			m_iteration_count = iter_cnt;
			if (print_output) output_data(func_name, root);
			return root; 
		}

		double Newton_Raphson::find_root(bool print_output)
		{
			std::string func_name = "NEWTON";
			check_all(func_name);
			double root = m_lft_brc;
			if( m_rht_brc != root ) //i.e. we are not using an initial guess
			{
				root = (fabs(m_function(m_lft_brc)) < fabs(m_function(m_rht_brc))) 
					? m_lft_brc : m_rht_brc;
			}
			double delta;
			size_t iter_cnt = 0;
			while(m_error > m_tolerance && iter_cnt++ < m_iter_max)
			{
				delta = -m_function(root)/m_derivative(root);
				root += delta;
				m_error = fabs(delta/root); 
			}
			//store iter_count for output
			m_iteration_count = iter_cnt;
			if (print_output) output_data(func_name, root);
			return root;
		}

		double Hybrid_B_N_R::find_root(bool print_output)
		{
			std::string func_name = "HYBRID-NEWTON";
			check_all(func_name);
			double root = m_rht_brc;
			double FA = m_function(m_lft_brc);
			double FR = m_function(root);
			double DFR = m_derivative(root);
			double delta;
			size_t iter_cnt = 0;
			while(m_error > m_tolerance && iter_cnt++ < m_iter_max)
			{
				if( (DFR*(root - m_lft_brc) - FR) * (DFR * (root - m_rht_brc) - FR) <= 0 )
				{
					delta = -FR/DFR;
					root += delta;
				}else
				{
					delta = (m_rht_brc - m_lft_brc)/2;
					root = (m_lft_brc + m_rht_brc)/2;
				}
				FR = m_function(root);
				DFR = m_derivative(root);
				if(FA*FR <= 0)
				{
					m_rht_brc = root;
				}else
				{
					m_lft_brc = root;
					FA = FR; 
				}
				m_error = fabs(delta/root);
			}
			//store iter_count for output
			m_iteration_count = iter_cnt;
			if (print_output) output_data(func_name, root);
			return root;
		}


		double Hybrid_B_S::find_root(bool print_output)
		{
			std::string func_name = "HYBRID-SECANT";
			check_all(func_name);
			double root = m_lft_brc;
			double delta;
			double FA = m_function(m_lft_brc);
			double FB = m_function(m_rht_brc);
			double FR = m_function(root);
			double DFR = (FB - FA)/(m_rht_brc - m_lft_brc);
			size_t iter_cnt = 0;
			while(m_error > m_tolerance && iter_cnt++ < m_iter_max)
			{
				if( (DFR*(root - m_lft_brc) - FR)*(DFR*(root - m_rht_brc) - FR) <= 0)
				{
					delta = -FR/DFR;
					root += delta;
				}else
				{
					delta = (m_rht_brc - m_lft_brc)/2;
					root = (m_lft_brc + m_rht_brc)/2;
				}

				FR = m_function(root);

				if(FA*FR <= 0)
				{
					m_rht_brc = root;
					FB = FR; 
				}else
				{
					m_lft_brc = root;
					FA = FR;
				}

				DFR = (FB - FA)/(m_rht_brc - m_lft_brc); 
				m_error = fabs(delta/root);
			}
			//store iter_count for output
			m_iteration_count = iter_cnt;
			if (print_output) output_data(func_name, root);
			return root;
		}


	}

}
