#ifndef ROOTSEARCH_HPP
#define ROOTSEARCH_HPP

#include <string>
#include <iostream>
#include <sstream>
#include <cstdio>

namespace phys {

/*
 *  Would it be better if the brackets for a root be passed to the find_root()
 *  interface function, rather than being initialised in the constructor, or
 *  explicitly set by the relevant setter functions? e.g.
 *
 *  	double find_root(double left, double right, bool print_output = false)
 *
 *  If so, how might this affect the internal representation of a RootSearch class?
 */

/** Virtual base class for the root search methods */
class RootSearch {
protected:
	/*	Default constructor */
	RootSearch();
	/*	Constructor using just the search m_function */
	RootSearch(double (*func)(double));
	/*	Typical constructor */
	RootSearch(double (*func)(double), double left, double right);
	/*	Specific constructor for the Newton-Raphson method with only the function and derivative */
	RootSearch(double (*func)(double), double (*deriv)(double));
	/* Specific constructor for the Newton-Raphson method */
	RootSearch(double (*func)(double), double (*deriv)(double),
			double initial_guess);
	/* Specific constructor for the hybrid Bisection-Newton-Raphson method */
	RootSearch(double (*func)(double), double (*deriv)(double), double left,
			double right);
public:
	/*	Destructor */
	virtual ~RootSearch() {
	}
private:
	/* No Copying */
	RootSearch(const RootSearch&);
	/* No Copying */
	const RootSearch& operator=(const RootSearch&);
public:
	//Interface function to find a root using a derived method
	virtual double find_root(bool print_output = false)=0;

	/*	Finds the brackets of a root between [start,end];
	 returns true on first bracket found, or false if no bracket found
	 start has to be less than end
	 */
	bool find_brackets(double start, double step, double end);

	void set_function(double (*func)(double)) {
		m_function = func;
	}

	void set_tolerance(double tol_to_set) {
		m_tolerance = tol_to_set;
	}

	//function pointer must be provided before using this function
	//- it checks a root is contained between the limits provided
	void set_brackets(double new_lft_brc, double new_rht_brc);

	//Newton-Raphson specific - use if you change the search function
	void set_derivative(double (*deriv)(double)) {
		m_derivative = deriv;
	}
	//Newton-Raphson specific - use if you want to set/change the initial guess for a root
	void set_initial_guess(double first_guess) {
		m_lft_brc = m_rht_brc = first_guess;
	}

	//Query the left bracket
	double lft_bracket() const {
		return m_lft_brc;
	}
	//Query the right bracket
	double rht_bracket() const {
		return m_rht_brc;
	}
	//Query the error achieved
	double get_error() const {
		return m_error;
	}
protected:
	/*	Checks that the m_function has been set, the brackets are set, and if a Newton method that the m_derivative has been set,
	 and/or if a bisection method that the brackets contain a root. */
	void check_all(std::string call_func);

	/* Outputs the result of a particular root search, including the error achieved, and the number of iteration*/
	void output_data(const std::string &calling_func, double root);
private:
	/*	Defaults: m_tolerance = 1.e-8; iteration_count = 0; iter_max = 100; error = 1.0; m_is_bisec = true;	m_is_newton = false; */
	void initialise_consts();
	void check_function_set(const std::string &calling_function);
	void check_derivative_set(const std::string &calling_function);
	void check_brackets_set(const std::string &call_func);
	void check_root_is_bracketed(const std::string &call_func, double func_lft,
			double func_rht);
protected:
	double (*m_function)(double);//!< pointer to m_function on which to root search.
	double (*m_derivative)(double);	//!< [optional] pointer to m_derivative of the m_function (Newton methods)
	double m_lft_brc;						//!< left hand bracket of a root
	double m_rht_brc;						//!< right hand bracket of a root
	double m_tolerance;					//!< Accuracy to achieve default 1.e-8
	double m_error;					//!< (Estimated) error in the root found.
	size_t m_iter_max;//!< Maximum number of iterations to try (prevents infinite loops) default 100
	size_t m_iteration_count;//!< Number of iterations performed to achieve m_tolerance.
	bool m_is_bisec;//!< Boolean flag; is this a bisection method? default true
	bool m_is_newton;//!< Boolean flag; is this a Newton method? default false
};

/** Bisection method */
class Bisection: public RootSearch {
public:
	Bisection() :
			RootSearch() {
	}
	Bisection(double (*func)(double)) :
			RootSearch(func) {
	}
	Bisection(double (*func)(double), double left, double right) :
			RootSearch(func, left, right) {
	}
	double find_root(bool print_output = false);
};

/** Secant method */
class Secant: public RootSearch {
public:
	Secant() :
			RootSearch() {
		m_is_bisec = false;
	}
	Secant(double (*func)(double)) :
			RootSearch(func) {
		m_is_bisec = false;
	}
	Secant(double (*func)(double), double first, double second) :
			RootSearch(func, first, second) {
		m_is_bisec = false;
	}
	double find_root(bool print_output = false);
};

/** Newton-Raphson method */
class Newton_Raphson: public RootSearch {
public:
	Newton_Raphson() :
			RootSearch() {
		m_is_bisec = false;
		m_is_newton = true;
	}
	Newton_Raphson(double (*func)(double), double (*deriv)(double)) :
			RootSearch(func, deriv) {
		m_is_bisec = false;
		m_is_newton = true;
	}
	Newton_Raphson(double (*func)(double), double (*deriv)(double),
			double first_guess) :
			RootSearch(func, deriv, first_guess) {
		m_is_bisec = false;
		m_is_newton = true;
	}
	double find_root(bool print_output = false);
};

/** Hybrid Bisection Newton-Raphson */
class Hybrid_B_N_R: public RootSearch {
public:
	Hybrid_B_N_R() :
			RootSearch() {
		m_is_newton = true;
	}
	Hybrid_B_N_R(double (*func)(double), double (*deriv)(double)) :
			RootSearch(func, deriv) {
		m_is_newton = true;
	}
	Hybrid_B_N_R(double (*func)(double), double (*deriv)(double), double left,
			double right) :
			RootSearch(func, deriv, left, right) {
		m_is_newton = true;
	}
	double find_root(bool print_output = false);
};

/** Hybrid Bisection Secant */
class Hybrid_B_S: public RootSearch {
public:
	Hybrid_B_S() :
			RootSearch() {
	}
	Hybrid_B_S(double (*func)(double)) :
			RootSearch(func) {
	}
	Hybrid_B_S(double (*func)(double), double left, double right) :
			RootSearch(func, left, right) {
	}
	double find_root(bool print_output = false);
};

} //namespace

#endif //header guard
