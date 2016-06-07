#ifndef ODE_SOLVERS_HPP
#define ODE_SOLVERS_HPP

#include "Storage.h"
#include "Differentials.h"
#include "DynVector.h"

#ifndef DBL_EPSILON
#define DBL_EPSILON 2.2204460492503131E-16
#endif

/*	This file contains the base class ODE_solver, and the derived solvers.
*/

namespace phys{
	namespace ode{

		class ODE_solver{
		protected:						
			/*	Constructor for the ODE_solver class. Protected as this is a virtual base class
				i.e. we cannot make a direct instance of it.
			*/
			ODE_solver(		phys::diffs::Diff_eqn* p_diff, 
							const state& system, 
							double step	);
		public:
			//virtual destructor as this is an abstract base class
			virtual ~ODE_solver(){}
		public:
			//interface
			/* Advances the system a single step.*/
			virtual state solve()=0;

			/*	Solves the system from initial independent value to end. */
			virtual phys::storage::ODE_Storage fullSolve(double end) = 0;

			/* Solves the system from initial independent value to end, wrapping the dependent variable to the range -pi to +pi*/
			virtual phys::storage::ODE_Storage fullSolveWrapped(double end)
			{
				size_t N = num_steps(m_current.x, end), count = 0;
				phys::storage::ODE_Storage container(m_current);
				while(count++ < N){
					this->wrap();
					container.store(m_current); //need to store the initial state before first solve
					this->solve();
				}
				wrap();
				container.store(m_current); //need to store the result of the final solve.
				return container;
			}

			/* Wraps the dependent values of the current step to the range -pi to +pi*/
			void wrap(){
				for (size_t i = 0; i < m_dims; i++){
					using namespace constants; //for PI
					if (m_current.y[i] < -PI) m_current.y[i] += 2 * PI;
					if (m_current.y[i] > +PI) m_current.y[i] -= 2 * PI;
				}
			}

			/*
				This function is primarily used to change the (initial) system state for a solver. 
				It can also be use to change the differential equation being solved, and/or the step size.

				To change the step size only use set_step(double). 

				There is no set_equation( Diff_eqn* )  function; if you're changing the equation it is very likely
				you'll want a different initial state. 
			*/
			void set_system(const state& system, phys::diffs::Diff_eqn* p_diff = 0, double step = double());


			/*	Use to change the step size	*/
			void set_step(double step);

			state get_current() const { return m_current;}
		protected:
			/*
				Prototype for solver methods that require an initialisation step.
				Derived methods that do not require initialisation call this function
				that is implicitly inlined and has a null body - it should be removed
				by compiler optimisation.
			*/
			virtual void init_method(){ /* null body */ }
			stdVec_d deriv(double x, const stdVec_d& y);
			stdVec_d deriv_B(double x, const stdVec_d& y = stdVec_d());
			size_t num_steps(double start, double& end);
		private:
			//function called in body of constructor
			void initialise_solver();

			//error message functions
			void __error_message_order_1() const;
			void __error_message_order_2() const;
		protected:
			//internal representation
			phys::diffs::Diff_eqn* m_ptr_diff_eqn;	//!< pointer to class that contains the equation to solve
			state m_current;						//!< system values x,y(,dy) updated at each step
			double m_step;							//!< step size for the integrator

			stdVec_d m_deriv_result;				//!< result of the derivative functions.
			stdVec_d m_k[6];						//!< typically used to store intermediate derivative values.

			size_t m_dims;							//!< dimensionality of the problem
			size_t m_num_of_vars;					//!< total number of variables in the problem
			size_t m_mid_idx;						//!< index at which y changes to dy in the dependent vector.
		};

		class Euler: public ODE_solver{
		public:
			Euler(		phys::diffs::Diff_eqn* p_diff, 
						const state& system, 
						double h	) 
				: ODE_solver(p_diff, system, h){}
			state solve();
			phys::storage::ODE_Storage fullSolve(double end);
		};

		class Imp_Euler: public ODE_solver{
		public:
			Imp_Euler(		phys::diffs::Diff_eqn* p_diff, 
							const state& system, 
							double h	) 
				: ODE_solver(p_diff, system, h){}
			state solve();
			phys::storage::ODE_Storage fullSolve(double end);
		};

		class Mod_Euler: public ODE_solver{
		public:
			Mod_Euler(		phys::diffs::Diff_eqn* p_diff, 
							const state& system, 
							double h	) 
				: ODE_solver(p_diff, system, h){}
			state solve();
			phys::storage::ODE_Storage fullSolve(double end);
		};

		class RK4: public ODE_solver{
		public:
			RK4(	phys::diffs::Diff_eqn* p_diff, 
					const state& system, 
					double h	) 
				: ODE_solver(p_diff, system, h){}
			state solve();
			phys::storage::ODE_Storage fullSolve(double end);
		};

		/*	On construction of a Stomer-Verlet solver the first step is automatically initialised.
			In order to save this data you need to set-up a (ODE) storage class object to .store() the 
			current system after construction and before the first .solve(). If using fullSolve
			all the data is stored automatically. Note the Stormer Verlet solver works only for differential
			equations of the type y'' = f(y), i.e. the differential function does not have an independent or
			first derivative term. 
		*/
		class Stormer_Verlet: public ODE_solver{
		public:
			Stormer_Verlet(	phys::diffs::Diff_eqn* p_diff, 
								const state& system, 
								double h	) 
				: ODE_solver(p_diff, system, h)
			{
				init_method();
			}
			state solve();
			phys::storage::ODE_Storage fullSolve(double end);
		private:
			virtual void init_method();
		private:
			state m_previous;
		};

		/*	Leapfrog integrates postion (dependent variable) at integer steps and the velocity 
			(first derivative) at half integer steps. First half step for velocity is computed
			using init_method() at construction (it uses an Euler step).
		*/
		class Leapfrog : public ODE_solver{
		public:
			Leapfrog(		phys::diffs::Diff_eqn* p_diff, 
							const state& system, 
							double h,
							bool synced = false) 
				:	ODE_solver(p_diff, system, h),
					m_synchronise(synced)
			{
				init_method();
			}
			state solve();
			phys::storage::ODE_Storage fullSolve(double end);
			
			/*	
				Use this function to synchronise the velocity to the position after each step 
				To store the advanced synced system in your main do: leapfrog->solve(); container.store(leapfrog->sync());
				The flag synchronise has no effect on this function.
			*/
			state sync();
			/*	
				Use this function after a full solve to synchronise the velcocity to the position 
				This can be selected to run automatically by switching the synchronise flag to true
			*/
			phys::storage::ODE_Storage sync ( const phys::storage::ODE_Storage& data );
		private:
			virtual void init_method();
		private:
			bool m_synchronise;	//!< is the velcoicty to be synchronised to the position or not? default false
		};

		class Velocity_Verlet: public ODE_solver{
		public:
			Velocity_Verlet(	phys::diffs::Diff_eqn* p_diff, 
								const state& system, 
								double h	) 
				: ODE_solver(p_diff, system, h){}
			state solve();
			phys::storage::ODE_Storage fullSolve(double end);
		};

		enum targetState { TARGET, NO_TARGET };

		class RKF45 : public ODE_solver{
		public:
			/*	Constructor for the adaptive step Runge-Kutta-Fehlberg solver.
				This will set up the required constants used in the method.
				Default values for error tolerance and hmin are as follows:
				1.e-5 and 1.e3 * DBL_EPSILON, respectively. (DBL_EPSILON ~ 2.e-16)
				hmax is set equal to the initial h you give in the constructor.
				Note that if going right to left then end < start, and h < 0.
			*/
			RKF45(	phys::diffs::Diff_eqn* p_diff, 
					const state& system, 
					double h_init,
					targetState target = NO_TARGET) : 
					ODE_solver(p_diff, system, h_init),
					alpha(9.e-1),
					tol(1.e-5),
					hmax(h_init),
					hmin(1.e3 * DBL_EPSILON),										
					htarg(hmax),
					m_targetState(target),
					m_targetHit(true),			//we take the initial state's x value as the first target
					m_nextTarget(system.x + htarg)
					{ init_constants();}
			void set_alpha(double alpha_to_set) {alpha = alpha_to_set;}
			void set_tolerance(double tol_to_set) {tol = tol_to_set;}
			void set_hmax(double hmax_to_set){hmax = hmax_to_set;}
			void set_hmin(double hmin_to_set){hmin = hmin_to_set;}			
			state solve();
			phys::storage::ODE_Storage fullSolve(double end);
			virtual phys::storage::ODE_Storage fullSolveWrapped(double end);

			bool isTargetHit(){return m_targetHit;}

		private:
			double alpha;
			double a1, a2, a3, a4, a5;
			double b10, b20, b21, b30, b31, b32, b40, b41, b42, b43, b50, b51, b52, b53, b54;
			double c0, c2, c3, c4, c5;
			double d0, d2, d3, d4, d5;
			double tol;
			double hmax, hmin, hnew, htarg;		
			double sign;	//!< indicates direction of itegration (l->r +ve, r->l -ve). RKF45 only.
			stdVec_d yhat;
			targetState m_targetState;
			bool m_targetHit;
			double m_nextTarget;
			void init_constants()
			{
				sign = (m_step < 0) ? -1.0 : 1.0;

				//am define the coefficients used to calculate x (the independent variable)
				a1 = 2.5e-1;
				a2 = 3.75e-1;
				a3 = 1.2e1 / 1.3e1;
				a4 = 1.e0;
				a5 = 5.e-1;

				//bmn define the coefficients used to calculate y (the dependent function)
				b10 = 2.5e-1;
				b20 = 3.e0 / 3.2e1;
				b21 = 9.e0 / 3.2e1;
				b30 = 1.932e3 / 2.197e3;
				b31 = -7.2e3 / 2.197e3;
				b32 = 7.296e3 / 2.197e3;
				b40 = 4.39e2 / 2.16e2;
				b41 = -8.e0;
				b42 = 3.68e3 / 5.13e2;
				b43 = -8.45e2 / 4.104e3;
				b50 = -8.e0 / 2.7e1;
				b51 = 2.e0;
				b52 = -3.544e3 / 2.565e3;
				b53 = 1.859e3 / 4.104e3;
				b54 = -1.1e1 / 4.e1;

				//cm are used to compute yhat
				c0 = 1.6e1 / 1.35e2;
				c2 = 6.656e3 / 1.2825e4;
				c3 = 2.8561e4 / 5.643e4;
				c4 = -9.e0 / 5.e1;
				c5 = 2.e0 / 5.5e1;

				//dm are used to compute the error in the integration step
				d0 = 1.e0 / 3.6e2;
				d2 = -1.28e2 / 4.275e3;
				d3 = -2.197e3 / 7.524e4;
				d4 = 1.e0 / 5.e1;
				d5 = 2.e0 / 5.5e1;
			}
		};

		/**	The Numerov class solves a single step of a 2nd order differential equation of the form
				u''(x) + q(x)u(x) = 0
			where we set our differential equation function to q(x). 
			To use, set up your system using the state constructor for a 2nd order ode with one dimension.
			See the example program that uses the Numerov method to solve the wavefunctions and energies 
			of a quantum particle bound in a harmonic potential (Harmonic_states.cpp).
		*/
		class Numerov : public ODE_solver
		{
		public:
			Numerov(	phys::diffs::Diff_eqn* p_diff, 
						const state& initial, 
						double h	):
						ODE_solver(p_diff, initial, h)
						{ init_method();}		
			state solve();
			phys::storage::ODE_Storage fullSolve(double end);
		private:
			virtual void init_method();
		private:
			void __error_order() const;
			void __error_dimensions() const;
		private:
			double m_h_sqr_12;
		};

		/**	The Numerov_s class solves a single step of a 2nd order differential equation of the form
				u''(x) + q(x)u(x) = s(x)
			where q(x) is our differential equation function and s(x) is the so called source function.
		*/
		class Numerov_s : public ODE_solver
		{
		public:
			Numerov_s(	phys::diffs::Diff_eqn* p_diff, 
							const state& initial, 
							double h, 
							double (*s)(double) ):
							ODE_solver(p_diff, initial, h), 
							m_src_func(s)	{ init_method();}
			void set_src_func(double (*s)(double)){ m_src_func = s; }			
			state solve();
			phys::storage::ODE_Storage fullSolve(double end);
		private:
			virtual void init_method();
		private:
			void __error_order() const;
			void __error_dimensions() const;
		private:
			double m_h_sqr_12;
			double(*m_src_func) (double);
			double m_s[3];
		};
	}
}

#endif 
