#ifndef DIFFERENTIALS_HPP
#define DIFFERENTIALS_HPP

#include "Maths.h"
#include "PhysicalUnits.h"

namespace phys{
	namespace diffs{

		class Diff_eqn{
		protected:
			/* Constructor */
			Diff_eqn(size_t eqn_order = 2) :order(eqn_order){ check_order(); }
		public:
			/* Destructor */
			virtual ~Diff_eqn(){}
		private:
			/* Prevent copying */
			Diff_eqn( const Diff_eqn& );
			const Diff_eqn& operator= (const Diff_eqn& );
		public:			
			/* Interface function to the ODE */
			virtual double differential_function(	double inde_var, 
													const stdVec_d& depend_vars,
													int num_of_dims, int var_flag)=0;
			void set_order(size_t order_to_set)
			{
				order = order_to_set;	
				check_order();
			}		
			int get_order() const {return order;}
		protected:
			size_t order;
		private:

			void check_order()
			{
				if (order > 2 ) 
				{
					std::string errmsg = "Only order one and order two ODEs supported ";
					throw(std::runtime_error(errmsg));
				}
			}
			
		};

		class SHM_eqn : public Diff_eqn{
		private:
			double k;						//!< Spring constant
			double D;						//!< Drag coefficient
			double(*drive_func)(double);	//!< Driving force function (user defined)
		public:
			/*	Default constructor
				Default values: spring const, k = 9.0; drag coeff, D = 0.0.
			*/
			SHM_eqn():	Diff_eqn(), 
						k(9.0), 
						D(0.0), 
						drive_func() {}
			/*	Constructor to set spring constant and drag.
			*/
			SHM_eqn(	double spring, 
						double drag		): 
						Diff_eqn(),
						k(spring), 
						D(drag), 
						drive_func(){}
			/*	Constructor for a user defined drive function */
			SHM_eqn(	double spring, 
						double drag, 
						double(*d_func)(double)		): 
						Diff_eqn(),
						k(spring), 
						D(drag),
						drive_func(d_func){}
			double differential_function(	double x, 
											const stdVec_d& y,
											int N, int c	)
			{
				double retval = -k*y[c] - D*y[c + N];
				if(drive_func != 0) retval += drive_func(x);
				return retval;
			}
			void set_spring_const(double new_k)	{k = new_k;}
			void set_drag_coeff(double new_D)	{D = new_D;}
			void set_drive_func(double(*f)(double)){drive_func = f;}
			void turn_off_drive(){drive_func = 0;}
		};

		class Gravity : public Diff_eqn{
		public:
			Gravity(	double drag = double(), 
						int turb = -1	) : 
						Diff_eqn(),
						D(drag), 
						turbulent(turb){}
			double differential_function(	double x, 
											const stdVec_d& y,
											int N, int c	) 
			{				
				double r;
				switch(N){
				case 1:
					r = fabs(y[0]);
					break;
				case 2:
					r = phys::maths::mag2D(y[0], y[1]);
					break;
				case 3:
					r = phys::maths::mag3D(y[0], y[1], y[2]);
					break;
				default:
					__error_gravity_dims();
				}
				using namespace constants; //for GS
				double result = -GS * y[c]/r/r/r;
				if( D != double())
				{
					double vel = y[c+N]; //sign of the velocity important
					if(turbulent) vel *= fabs(y[c+N]); //use drag proportional to v squared.
					result -= D*vel; //minus as drag force acts in opposite direction to velocity.
				}
				return result;
			}
			void set_drag_coeff(double drag){D=drag;}
			void toggle_turbulence(){turbulent = -turbulent;}
		private:
			double D;					//!< set to non-zero value to add air resistance.
			int turbulent;				//!< set to true to use v squared model of drag.

			void __error_gravity_dims()
			{
				std::string errmsg = "Gravity class only supports 1, 2, or 3 spatial dimensions";
				throw (std::runtime_error(errmsg));
			}
		};

		class Van_Der_Pol : public Diff_eqn{
		public:
			Van_Der_Pol() : Diff_eqn(), 
							mu(1.0) {}
			Van_Der_Pol(double mu_val) :	Diff_eqn(), 
											mu(mu_val) {}
			double differential_function(double x, const stdVec_d& y, int N, int c)
			{
				return mu*(1 - y[c]*y[c])*y[c+N] - y[c];
			}
			void set_mu(double new_mu){mu = new_mu;}
		private:
			double mu;			//!< Damping parameter. Default = 1.0
		};

		class Duffing : public Diff_eqn{
		public:
			/*	Default constructor sets up equation as SHM i.e. only beta non zero (= 9.0) */
			Duffing() : Diff_eqn(), 
						beta(9.0), 
						alpha(0.0), 
						delta(0.0), 
						gamma(0.0), 
						omega(0.0){}
			//----------------------------------------------------------------------------------
			Duffing(	double stiffness, 
						double nonLinear, 
						double damping, 
						double amplitude, 
						double frequency	) :
						Diff_eqn(),
						beta(stiffness), 
						alpha(nonLinear), 
						delta(damping), 
						gamma(amplitude), 
						omega(frequency){}
			double differential_function(	double x, 
											const stdVec_d& y,
											int N, int c	)
			{
				return (gamma*cos(omega*x)-(beta*y[c] + alpha*y[c]*y[c]*y[c] + delta*y[c+N])); 
			}
			void set_alpha(double a){alpha = a;} //!< non-linear
			void set_beta (double b){beta  = b;} //!< stiffness
			void set_gamma(double c){gamma = c;} //!< drive amplitude
			void set_delta(double d){delta = d;} //!< damping
			void set_omega(double w){omega = w;} //!< drive frequency
		private:
			double beta;		//!< Stiffness (spring const.)
			double alpha;		//!< Non-linear restoring force coefficient			
			double delta;		//!< Damping coefficient
			double gamma;		//!< Amplitude of the periodic driving force			
			double omega;		//!< Frequency of the driving force
		};


		/*	Use to define your own differential_function using the member prototype. Default
			order is 2. Use set_order to change to 1 if necessary. Note that only order
			1 and 2 ODEs supported. The vector 'params' is available for you to use in 
			your differential_function; just use set_parameters function in your main.
		*/
		class User_eqn : public Diff_eqn{
		public:
			void set_parameters(const stdVec_d& user_parameters)
			{
				params = user_parameters;
			}
			/*	@brief All arguments automatically passed in by the ODE_solver of choice.
				@param x double value of the independent varaible
				@param y vector<double> of the state of the system (dependent values in first half, first derivative in the second half)
				@param N integer variable of the number of dimensions in the system
				@param c integer variable of the specific dimension for which to compute the differential
			*/
			double differential_function(	double independent, 
											const stdVec_d& dependent,
											int num_of_dims, int dim_choice);
		private:
			stdVec_d params;
		};

	}
}
#endif
