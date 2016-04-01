#ifndef QUADRATURE_HPP
#define QUADRATURE_HPP

#include <string>
#include <iostream>
#include <sstream>
#include <cstdio>
#include <exception>

#include "DynVector.h"


namespace phys{
	namespace quad{

		/*	Virtual base class for a non-Gaussian quadrature object */
		class Quadrature{	
		protected:
			//Construtor to specify the accuracy order
			Quadrature(size_t n = 2);
		public:
			virtual ~Quadrature(){}
		private:
			/* No Copy */
			Quadrature(const Quadrature&);
			/* No Copy */
			const Quadrature& operator= (const Quadrature&);
		public:
			//Interface function to perform the derived method of integration:
			//@param	*func pointer to the integrand function (pass the function identifier as the argument)
			//			note integrand takes one double argument only.
			//@param	lft double value of the left hand limit of the integration
			//@param	rht double value of the right hand limit of the integration
			//@param	strips unisgned int of the number of strips to use in the composite quadrature method
			//@return	the value of the numerical quadrature for the parameters specified.
			virtual double integrate(double(*func)(double), double lft, double rht, size_t strips) = 0; 

			//Query functions
			double get_error() const {return m_error;}
			size_t get_order() const {return m_order;}
			size_t get_func_calls() const { return m_func_calls; }

			void reset_error() { m_error = 1.; }
			void reset_func_calls(){ m_func_calls = 0; }

			//interface function to return the name of the method used
			virtual std::string get_name()=0;
		protected:
			double m_error;		//!< estimated error in the numerical quadrature.
			size_t m_order;		//!< order of the error improvement for a particular method (e.g. Simpson is order 4)
			size_t m_func_calls;	//!< number of function calls made during integration (use for assessment of different methods)
		};

		/* Mid-Ordinate method. Derived class of Quadrature*/
		class Mid_ordinate : public Quadrature{
		public:
			Mid_ordinate() : Quadrature(){}
			double integrate(double(*f)(double), double lft, double rht, size_t strips);
			std::string get_name(){return "Mid";}
		};

		/* Trapeziod method. Derived class of Quadrature*/
		class Trapeziod : public Quadrature{
		public:
			Trapeziod() : Quadrature() {}
			double integrate(double(*f)(double), double lft, double rht, size_t strips);
			std::string get_name(){return "Trap";}
		};
		
		/* Simpson method. Derived class of Quadrature*/
		class Simpson : public Quadrature{
		public:
			Simpson() : Quadrature(4) {}
			double integrate(double(*f)(double), double lft, double rht, size_t strips);
			std::string get_name(){return "Simp";}
		};
		
		/*	Boole's method. Derived class of Quadrature.
			To estimate the error in Boole we would have to find finite difference scheme for the 6th 
			derivative. Instead, we estimate the error by first computing the integral for a step size of 
			h, s_h. Then, doubling the number of strips we compute s_2h. Using the knowledge that Boole's 
			rule is 6th order accurate, we can then compute the error as
			err_2h = |s_h - s_2h|/63
			The returned integral is then s_2h +/- err_2h. You select this error computation by setting 
			error_wanted to true; it is false by default and, as such, Boole will only return s_h with no 
			estimate of the error.  
		*/
		class Boole : public Quadrature{
		public:
			Boole() :Quadrature(6), error_wanted(false) {}
			Boole(bool err_req) :Quadrature(6), error_wanted(err_req) {}
			double integrate(double(*f)(double), double lft, double rht, size_t strips);
			std::string get_name(){return "Boole";}
		private:
			bool error_wanted;		//!< toggle for error estimation at 2*strips.
		};

		/********************************************************************************************* 
		*	Gaussian Quadrature 
		**********************************************************************************************/
		class Gauss{
		protected:
			//Constructor: specify the number of knots to use in the Gaussian quadrature
			Gauss(size_t pts);
		public:
			virtual ~Gauss(){}
		private:
			/* No Copy */
			Gauss( const Gauss& );
			/* No Copy */
			const Gauss& operator=( const Gauss& ); 
		public:
			//Interface function to perform the derived Gaussian method of integration
			//@param	*func pointer to the integrand function (pass the function identifier as the argument)
			//			note integrand takes one double argument only.
			//@param	lft double value of the left hand limit of the integration
			//@param	rht double value of the right hand limit of the integration (for semi-infinite pass an empty double)
			//@return	the numerical solution of the quadrature.
			virtual double integrate(double(*func)(double), double lft = 0., double rht = double() ) = 0;
			/*
			Use this function to change the number of points to use in this instance of the quadrature.
			*/
			void set_points(size_t pts);
			/*
			Use this function to change the precision of COMPUTED abscissas and weights.
			Only utilised when the number of points wanted is not on the predefined list.
			Predefined points are: 2 - 20, 32, 64, 96, 100, 128, 256, 512, 1024.
			*/
			void set_precision(double acc);

			size_t get_func_calls() const { return m_func_calls; }
			void reset_func_calls(){ m_func_calls = 0; }
		protected:
			//This populates the knots and corresponding weights to use in the quadrature
			//called in the body of the Gauss constructor, or when you want to change the
			//number of points, or the precision of the COMPUTED knots and weights. 
			virtual void initialise() = 0;
			//This function COMPUTES the knots and weights for the number of points speficied
			//in the quadrature, iff these values are not predefined in "Gauss_knots_weights.hpp"
			virtual std::pair<stdVec_d, stdVec_d> compute_x_w() = 0;
		protected:
			size_t m_numKnots;				//!< Number of knots to use in the quadrature
			stdVec_d m_knots;				//!< Vector of the knots to use in the quadrature
			stdVec_d m_weights;				//!< Vector of the corresponding weights.
			size_t m_func_calls;			//!< Number of function calls made to compute solution
			double m_precision;				//!< Accuracy for the computed knots and weights; Default is 1.e-10										
		};

	
		/*	Class for the Gauss-Legendre quadrature. This is a derived class of Gauss.
			Note that this will give exact solutions to polynomial functions of order 2n-1 or less.
		*/
		class Legendre: public Gauss {
		public:
			/*	
				Constructor for the Legendre class. Use this if you intend to keep the points
				in the quadrature constant. The initialise function will either fetch the required
				knots and weights from a predefined list, or will compute them if not defined.
				Predefined points: 2-20,16,32,64,96,100,128,256,512,1024.
			*/
			Legendre(size_t points = 2);
			/*	
				This function performs the integration.
			*/
			virtual double integrate(double(*f)(double), double lft, double rht );
			std::string get_name(){return "Legendre";}
		private:
			virtual void initialise();
			virtual std::pair< stdVec_d, stdVec_d >compute_x_w();
		};	
	
		class Laguerre : public Gauss {
		public:
			/*	
				Constructor for the Laguerre class. If the weight function is implict in your integrand
				set modify to true;this will multiply the weights by the exponential of the 
				corresponding abscissa. To explain:

				If your intgral is
				|a,+inf) exp(-x)*f(x)dx, then set modify false; else
				|a,+inf) f(x)dx  then set modify true (this is default).

				where f(x) is the function you pass to the integrate call; 'a' can be non-zero.
				To check the correct set up note that |0,+inf) exp(-x) * x*x*x == 6. 

			*/
			Laguerre(size_t points = 1, bool modify=true);		
			/*	
				This function performs the integration. The left limit is zero by default. 
			*/
			virtual double integrate(double(*f)(double), double lft = 0., double = double() );
			std::string get_name(){return "Laguerre";}
		private:
			bool modified;						//!< is the exponetial implicit or not?
			virtual void initialise();
			std::pair<stdVec_d, stdVec_d> compute_x_w();
		};

		/*******************************************************************************************************
		*	Other methods
		*******************************************************************************************************/

		/*	Adaptive class. This is NOT a derived class of Quadrature however we can pass it a pointer to a
			derived class of quadrature for adapting; the default is the trapezoid method. The class implements 
			a recursive function that subsequently halves the size of the integration interval until an error 
			tolerance is met.
		*/
		class Adaptive {
		public:
			Adaptive() :	m_pQuad(new Trapeziod), 
							m_I(), 
							m_order(m_pQuad->get_order()),
							m_tolerance(1.e-6), 
							m_count(0), 
							m_cnt_min(1), 
							m_cnt_max(100),
							m_num_f_calls(0),
							m_print_output(false),
							d_ctor(true){}

			Adaptive(	Quadrature* ptr_q, 
						double tol = 1.e-6,
						bool print_o = false	) : 
						m_pQuad(ptr_q), 
						m_I(), 
						m_order(m_pQuad->get_order()),
						m_tolerance(tol), 
						m_count(0), 
						m_cnt_min(1), 
						m_cnt_max(100),
						m_num_f_calls(0),
						m_print_output(print_o),
						d_ctor(false) {}

			~Adaptive()
			{
				if(d_ctor)
					delete m_pQuad;
			}

			double integrate(double(*f)(double), double lft, double rht, size_t strips = 1, double tol = double());
			void reset_count(){m_count = 0;}
			void set_quadrature_method(Quadrature* ptr_q)
			{
				m_pQuad = ptr_q;
				m_order = m_pQuad->get_order();
				reset_count();
			}
			void set_tolerance(double tol)
			{
				m_tolerance = tol;
				reset_count();
			}
			void toggle_output(){m_print_output = (m_print_output == false) ? true:false;}
			int get_count() const {return m_count;}
			int get_f_calls() const{return m_num_f_calls;}

		private:
			double recursive(double lft, double rht, size_t strips, double tol, size_t&cnt, double&error);

		private:
			Quadrature* m_pQuad;			//!< pointer to the quadrature method to use
			double(*m_I)(double);				//!< integrand function (assigned in integrate function)
			size_t m_order;					//!< accuracy order (obtained from quad method). For Gaussian types == 1.
			double m_tolerance;				//!< precision to achieve.
			size_t m_count;					//!< No of times the adaptive method has halved the interval
			size_t m_cnt_min;					//!< Min no of times to halve to avoid any unexpected initial behaviour
			size_t m_cnt_max;					//!< Max no of times to halve the interval to avoid infinite recursions.
			size_t m_num_f_calls;				//!< Number of function calls made during integration.
			bool m_print_output;				//!< option to print output for each integration.
			bool d_ctor;					//!< flag to identify when the default constructor has been used.			
		};

		/*	Romberg's method of quadrature. This is NOT a derived class of Quadrature.
			It applies Richardson Extrapolation to the Trapeziod method to obtain a 
			numerical result to a desired accuracy; default 1.e-6.
		*/
		class Romberg {
		public:
			Romberg() :		m_tol(1.e-6), 
							m_level(10), 
							m_error(1.0), 
							m_print_output(false),
							m_func_calls(0){}
			Romberg( double acc, size_t lvl=10, bool print_o = false) : 
							m_tol(acc), 
							m_level(lvl), 
							m_error(1.0), 
							m_print_output(print_o),
							m_func_calls(0){}

			double integrate(double(*f)(double), double lft, double rht,
				double tolerance=double(), size_t levels=size_t());
 
			void set_accuracy(double accuracy) {m_tol = accuracy;}
			void set_max_level(size_t max_lvl){ m_level = max_lvl; }
			void toggle_output(){m_print_output = (m_print_output == false) ? true:false;}
			double get_error() const {return m_error;}
			size_t get_func_calls() const {return m_func_calls;}
		private:
			double m_tol;						//!< Required accuracy (upper bound). default = 1.e-6
			size_t m_level;					//!< Max extrapolation level to reach. default = 10
			double m_error;					//!< Actual error achieved (difference between the two most accurate extrapolations).
			bool m_print_output;				//!< Option to print output for each integration.
			size_t m_func_calls;				//!< Number of function calls used to get solution.
		}; 
	}
}


#endif
