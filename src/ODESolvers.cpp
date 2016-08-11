#include <cmath>

#include "ODESolvers.h"


/*	.solve() performs a single integration step of the solver used, and returns a state.
	It is used if we wish to terminate the integration on some other criteria rather than
	integrating over a fixed range of the independent variable. (e.g. measuring when a 
	damped oscillator has an amplitude less than 1% of it's starting value, say). Each solve
	function has the following layout:

	*::solve(){
		state next = m_current;
		//algorithm loop
		.
		.
		//update m_current at end of loop
		m_current = next
		return next;
	}

	.fullSolve() performs a full integration from initial independent value to the end
	value specified as an argument. It returns a storage object that contains all the steps 
	computed.The full solve takes a double value argument that is the end point of the 
	integration. Each fullSolve function (with the execption of the adaptive RKF, explained 
	shortly) then checks to see that the value (end - start)/step is an integer value 
	(it uses this value for condition checking in the integration loop). If this isn't an 
	integer value it adjusts the end point such that the value is the next largest integer, and 
	writes a message to the console that the end point has been adjusted and to what value. 
	We go to the next largest integer so that the end point the user specified is 
	contained in the integration.
		
	For the adaptive Runge-Kutta this integer value isn't defined by its very nature. For
	the loop condition checking we have to use the double end point value, which may mean 
	we go beyond the end point value but is this a problem? ... possibly if we need to hit
	a particular target. 

	With performance in mind fullSolve does not call solve() but implements the same algorithm
	within a while loop explicitly (removes unnecessary function calls), as well as storing the 
	data. fullSolve has the following structure for each method (apart from the adaptive method
	which doesn't perform the check at the start). 
	*::fullSolve(double end)
	{
		check the end value is an integer number of steps from the start value; adjust otherwise.
		initialise a storage object for our data.
		assign state m_current to state next (this initialises a tempory/local state holder).
		while loop to execute the algorithm till end - stores data in storage object
		stores final state to storage object outside loop (loop exits before we store the final state)
		return the storage object.
	}

	Note that because we have defined the necessary operators in the header physVector.hpp we can
	replace the for loops that scan over the entire vector m_current.y with just the operator syntax,
	removing the subscript i. To illustrate in the RK4 alogrithm we have the loop:

	for(size_t i = 0; i < m_num_of_vars; ++i)
		next.y[i] = m_current.y[i] + (m_step/2.)*m_k[0][i];

	which can be replaced with:

	next.y = m_current.y + (m_step/2.) * m_k[0];

	which reads better but is it more computationally efficient? 
*/

namespace phys{
	namespace ode{
		/* ---------------- ODE Virtual Base Class implementation ---------------------- */
		/*	Constructor for the ODE virtual base class. 
			Any member not in list is assigned using initialise_solver(), 
			except m_k which is assigned in the implementation of solve or fullSolve*/
		ODE_solver::ODE_solver(	phys::diffs::Diff_eqn* p_diff,
								const state& system,
								double step) :								
								m_ptr_diff_eqn(p_diff),
								m_current(system),
								m_step(step)
		{
			initialise_solver();
		}

		void ODE_solver::set_system(const state& system, phys::diffs::Diff_eqn* p_diff, double step)
		{
			m_current = system;
			if (p_diff != 0)
				m_ptr_diff_eqn = p_diff;
			if (step != double())
				m_step = step;
			initialise_solver();
			init_method();
		}

		void ODE_solver::set_step(double step)
		{ 
			m_step = step; 
			init_method(); 
		}

		void ODE_solver::initialise_solver()
		{
			int order = m_ptr_diff_eqn->get_order();
			if (order == 1 && m_current.dy.empty() == false)
				__error_message_order_1();
			if (order == 2 && m_current.dy.empty() == true)
				__error_message_order_2();
			m_dims = m_current.num_dims;
			m_num_of_vars = m_dims*order;
			m_deriv_result = stdVec_d(m_num_of_vars);
			m_mid_idx = (order == 1) ? 0 : m_dims;
		}

		/*	We split the deriv function into two so that we can autonomously handle both 1st and 
			2nd ordered ODEs without conditional branching statements. If 1st order them m_mid_idx = 0 
			and we skip the loop in deriv and only call deriv_B. Else if 2nd order m_mid_idx = num of 
			dimensions in the problem and necessarily shuffle the values in the m_deriv_result vector.
		*/
		stdVec_d ODE_solver::deriv(double x, const stdVec_d& y)
		{
			for(size_t i = 0; i < m_mid_idx; ++i)
				m_deriv_result[i]= y[i + m_mid_idx];			
			return deriv_B(x,y); 
		}

		stdVec_d ODE_solver::deriv_B(double x, const stdVec_d& y)
		{
			for(int i = m_mid_idx; i < m_num_of_vars; ++i)
				m_deriv_result[i] = m_ptr_diff_eqn->differential_function(x,y,m_dims,i-m_mid_idx);
			return m_deriv_result;
		}

		//@query: The tolerance for an end being close enough is hard coded - should this be changed to a variable?
		size_t ODE_solver::num_steps(double start, double& end)
		{			
			double N = (end - start)/m_step; 
			assert(N > 0); //Could we handle this better?

			double tN = floor(N); 

			if( N - tN < 1.e-6 ) return size_t(N); //the number of steps is close enough to an int value
			else{
				size_t retval = size_t(tN + 1.); //next largest integer
				end = retval*m_step + start; //adjust end point to reflect change.							
				return retval;
			}
		}

		/* --- Euler --- */
		/* Here we can directly update m_current w/o a local state being instantiated */
		state Euler::solve() 
		{
			m_k[0] = deriv(m_current.x, m_current.y); 
			for(size_t i = 0; i < m_num_of_vars; ++i)
				m_current.y[i] += m_step*m_k[0][i];
			m_current.x += m_step;
			return m_current;
		}

		phys::storage::ODE_Storage Euler::fullSolve(double end)
		{
			size_t N = num_steps(m_current.x, end), count = 0;
			phys::storage::ODE_Storage container(m_current);
			while( count++ < N )
			{
				container.store(m_current);
				m_k[0] = deriv(m_current.x, m_current.y); 			
				for(size_t i = 0; i < m_num_of_vars; ++i)
					m_current.y[i] += m_step*m_k[0][i];
				m_current.x += m_step;
			}
			container.store(m_current);
			return container;
		}

		/* --- Improved Euler --- */
		state Imp_Euler::solve()
		{
			state next = m_current;
			m_k[0] = deriv(next.x, next.y);
			for(size_t i = 0; i < m_num_of_vars; ++i)
				next.y[i] = m_current.y[i] + (m_step/2.)*m_k[0][i];
			m_k[1] = deriv(next.x + m_step/2., next.y);
			for(size_t i = 0; i < m_num_of_vars; ++i)
				next.y[i] = m_current.y[i] + m_step*m_k[1][i];
			next.x = m_current.x + m_step;
			m_current = next;
			return next;
		}

		phys::storage::ODE_Storage Imp_Euler::fullSolve(double end)
		{
			size_t N = num_steps(m_current.x, end), count = 0;
			phys::storage::ODE_Storage container(m_current);
			state next = m_current;
			while( count++ < N )
			{
				container.store(m_current);
				m_k[0] = deriv(next.x, next.y);
				for(size_t i = 0; i < m_num_of_vars; ++i)
					next.y[i] = m_current.y[i] + (m_step/2.)*m_k[0][i];
				m_k[1] = deriv(next.x + m_step/2., next.y);
				for(size_t i = 0; i < m_num_of_vars; ++i)
					next.y[i] = m_current.y[i] + m_step*m_k[1][i];
				next.x = m_current.x + m_step;
				m_current = next;
			}
			container.store(m_current);
			return container;
		}

		/* --- Modified Euler --- */
		state Mod_Euler::solve()
		{
			state next = m_current; 
			m_k[0] = deriv(next.x, next.y);
			for(size_t i = 0; i < m_num_of_vars; ++i)
				next.y[i] = m_current.y[i] + m_step*m_k[0][i];
			m_k[1] = deriv(next.x + m_step, next.y);
			for(size_t i = 0; i < m_num_of_vars; ++i)
				next.y[i] = m_current.y[i] + (m_step/2.)*(m_k[0][i] + m_k[1][i]); 
			next.x = m_current.x + m_step;
			m_current = next;
			return next;
		}

		phys::storage::ODE_Storage Mod_Euler::fullSolve(double end)
		{
			size_t N = num_steps(m_current.x, end), count = 0;
			phys::storage::ODE_Storage container(m_current);
			state next = m_current;
			while( count++ < N )
			{
				container.store(m_current);
				m_k[0] = deriv(next.x, next.y);
				for(size_t i = 0; i < m_num_of_vars; ++i)
					next.y[i] = m_current.y[i] + m_step*m_k[0][i];
				m_k[1] = deriv(next.x + m_step, next.y);
				for(size_t i = 0; i < m_num_of_vars; ++i)
					next.y[i] = m_current.y[i] + (m_step/2.)*(m_k[0][i] + m_k[1][i]); 
				next.x = m_current.x + m_step;
				m_current = next;
			}
			container.store(m_current);
			return container;
		}


		/* --- Runge-Kutta --- */
		state RK4::solve()
		{
			state next = m_current;
			//compute intermediate derivatives
			m_k[0] = deriv(next.x, next.y); 
			for(size_t i = 0; i < m_num_of_vars; ++i)
				next.y[i] = m_current.y[i] + (m_step/2.)*m_k[0][i];
			m_k[1] = deriv(next.x + m_step/2., next.y);
			for(size_t i = 0; i < m_num_of_vars; ++i)
				next.y[i] = m_current.y[i] + (m_step/2.)*m_k[1][i];
			m_k[2] = deriv(next.x + m_step/2., next.y);
			for(size_t i = 0; i < m_num_of_vars; ++i)
				next.y[i] = m_current.y[i] + m_step*m_k[2][i];
			m_k[3] = deriv(next.x + m_step, next.y);
			//Advance one step
			for(size_t i = 0; i < m_num_of_vars; ++i)
				next.y[i] = m_current.y[i] + (m_step/6.)*(m_k[0][i] + 
					2.*m_k[1][i] + 2.*m_k[2][i] + m_k[3][i]);
			next.x = m_current.x + m_step;
			m_current = next;
			return next;
		}

		phys::storage::ODE_Storage RK4::fullSolve(double end)
		{
			size_t N = num_steps(m_current.x, end), count = 0;
			phys::storage::ODE_Storage container(m_current);
			state next = m_current;
			while( count++ < N )
			{
				container.store(m_current);
				m_k[0] = deriv(next.x, next.y);
				for (size_t i = 0; i < m_num_of_vars; ++i)
					next.y[i] = m_current.y[i] + (m_step/2.)*m_k[0][i];
				m_k[1] = deriv(next.x + m_step/2., next.y);
				for (size_t i = 0; i < m_num_of_vars; ++i)
					next.y[i] = m_current.y[i] + (m_step/2.)*m_k[1][i];
				m_k[2] = deriv(next.x + m_step/2., next.y);
				for (size_t i = 0; i < m_num_of_vars; ++i)
					next.y[i] = m_current.y[i] + m_step*m_k[2][i];
				m_k[3] = deriv(next.x + m_step, next.y);
				for (size_t i = 0; i < m_num_of_vars; ++i)
					next.y[i] = m_current.y[i] + (m_step/6.)*(m_k[0][i] +
						2.*m_k[1][i] + 2.*m_k[2][i] + m_k[3][i]);
				next.x += m_step;
				m_current = next;
			}
			container.store(m_current);
			return container;
		}

		/* --- Stormer-Verlet --- */
		void Stormer_Verlet::init_method()
		{
			//need to compute first step in-order to advance the integrator.
			//assign the initialising system to m_previous
			m_previous = m_current;
			//compute first step from initial conditions
			m_k[0] = deriv(m_current.x, m_current.y);
			for (size_t i = 0; i < m_dims; i++)
				m_current.y[i] += m_step*m_k[0][i] + (m_step*m_step/2.)*m_k[0][i + m_dims];
			m_current.x += m_step;
		}

		//velocity is not updated with the Stormer-Verlet integrator.
		state Stormer_Verlet::solve()
		{
			state next = m_current;
			m_k[0] = deriv_B(m_current.x, m_current.y); 
			for(size_t i=0; i<m_dims; i++)
				next.y[i] = 2.*m_current.y[i] - m_previous.y[i] + (m_step*m_step)*m_k[0][i+m_dims];
			next.x = m_current.x + m_step;
			//assign m_current state to m_previous state for the next step.
			m_previous = m_current;
			m_current = next;
			return next;
		}

		phys::storage::ODE_Storage Stormer_Verlet::fullSolve(double end)
		{
			size_t N = num_steps(m_current.x, end), count = 0;
			phys::storage::ODE_Storage container(m_current);
			container.store(m_previous); //initial system data
			state next = m_current;
			while( count++ < N )
			{
				container.store(m_current);
				m_k[0] = deriv_B(m_current.x, m_current.y); 
				for(size_t i=0; i<m_dims; i++)
					next.y[i] = 2.*m_current.y[i] - m_previous.y[i] + (m_step*m_step)*m_k[0][i+m_dims];
				next.x = m_current.x + m_step;
				m_previous = m_current;
				m_current = next;
			}
			container.store(m_current);
			return container;
		}

		/* --- Leapfrog --- */ 
		void Leapfrog::init_method()
		{
			//compute acceleration from initial position
			m_k[0] = deriv_B(m_current.x, m_current.y);
			//compute the first half integer step for velocity
			for (size_t i = 0; i < m_dims; i++)
				m_current.y[i + m_dims] += (m_step/2.)*m_k[0][i + m_dims];
		}

		state Leapfrog::solve()
		{
			state next = m_current;
			for(size_t i = 0; i < m_dims; i++)
				next.y[i] = m_current.y[i] + m_step*m_current.y[i+m_dims];
			m_k[0] = deriv_B(next.x, next.y);
			for(size_t i = 0; i < m_dims; i++)
				next.y[i+m_dims] = m_current.y[i+m_dims] + m_step*m_k[0][i+m_dims];
			next.x = m_current.x + m_step;
			m_current = next;
			return next;
		}

		state Leapfrog::sync()
		{
			state synced_system = m_current;
			m_k[0] = deriv_B(m_current.x, m_current.y);
			for (size_t i = 0; i < m_dims; ++i)
				synced_system.y[i + m_dims] -= (m_step/2.)*m_k[0][i + m_dims];
			return synced_system;
		}

		/*	This will store the data with the velocity synchronised to the position as an option. */
		phys::storage::ODE_Storage Leapfrog::fullSolve(double end)
		{
			size_t N = num_steps(m_current.x, end), count = 0;
			phys::storage::ODE_Storage container(m_current);
			state next = m_current;
			while( count++ < N )
			{
				container.store(m_current);
				for(size_t i = 0; i < m_dims; i++)
					next.y[i] = m_current.y[i] + m_step*m_current.y[i+m_dims];
				m_k[0] = deriv_B(next.x, next.y);
				for(size_t i = 0; i < m_dims; i++)
					next.y[i+m_dims] = m_current.y[i+m_dims] + m_step*m_k[0][i+m_dims];
				next.x = m_current.x + m_step;
				m_current = next;
			}
			container.store(m_current);
			return (m_synchronise) ? sync(container) : container;
		}

		phys::storage::ODE_Storage Leapfrog::sync(const phys::storage::ODE_Storage& data)
		{
			//grab the stored system data
			stdVec_d inde = data.get_independent();
			stdVec_d depn = data.get_dependent();
			stdVec_d deri = data.get_first_deriv();

			phys::storage::ODE_Storage retval(state(inde[0], depn[0], deri[0]));
			size_t n = inde.size();
			//sync the velocity at each integration step
			for (size_t p = 0; p < n; ++p)
			{
				m_current = state(inde[p], depn[p], deri[p]);
				//store the synched system.
				retval.store(sync());
			}
			return retval;
		}


		/*--- Velocity-Verlet --- */
		state Velocity_Verlet::solve()
		{
			state next = m_current;
			//advance positions
			m_k[0] = deriv(next.x, next.y);		
			for(size_t i = 0; i < m_dims; ++i)
				next.y[i] = m_current.y[i] + m_step*m_k[0][i] + (m_step*m_step/2.)*m_k[0][i+m_dims];
			//advance velocity
			m_k[1] = deriv_B(next.x, next.y);			
			for(size_t i = 0; i < m_dims; ++i)
				next.y[i+m_dims] = m_current.y[i+m_dims] + (m_step/2.)*(m_k[0][i+m_dims] + m_k[1][i+m_dims]);
			next.x = m_current.x + m_step;
			m_current = next;
			return next;
		}

		phys::storage::ODE_Storage Velocity_Verlet::fullSolve(double end)
		{
			size_t N = num_steps(m_current.x, end), count = 0;
			phys::storage::ODE_Storage container(m_current);
			state next = m_current;
			while( count++ < N )
			{
				container.store(m_current);
				m_k[0] = deriv(next.x, next.y);		
				for(size_t i = 0; i < m_dims; ++i)
					next.y[i] = m_current.y[i] + m_step*m_k[0][i] + (m_step*m_step/2.)*m_k[0][i+m_dims];
				m_k[1] = deriv_B(next.x, next.y);			
				for(size_t i = 0; i < m_dims; ++i)
					next.y[i+m_dims] = m_current.y[i+m_dims] + (m_step/2.)*(m_k[0][i+m_dims] + m_k[1][i+m_dims]);
				next.x = m_current.x + m_step;
				m_current = next;
			}
			container.store(m_current);
			return container;			
		}

		/* --- Runge-Kutta-Fehlberg --- */
		state RKF45::solve()
		{
			//the initial state should be stored before calling this function

			state next = m_current;
			//assign m_current y vector to yhat; values will be overwritten later on.
			yhat = m_current.y;
			
			m_targetHit = false;

			while(next.x == m_current.x)//i.e. integrator not advanced
			{
				//compute intermediate values for derivatives
				m_k[0] = deriv(next.x, next.y);
				for(size_t i = 0; i < m_num_of_vars; ++i)
					next.y[i] = m_current.y[i] + b10*m_step*m_k[0][i];

				m_k[1] = deriv(next.x + a1*m_step, next.y);
				for(size_t i = 0; i < m_num_of_vars; ++i)
					next.y[i] = m_current.y[i] + m_step*(b20 * m_k[0][i] + b21*m_k[1][i]);

				m_k[2] = deriv(next.x + a2*m_step, next.y);
				for(size_t i = 0; i < m_num_of_vars; ++i)
					next.y[i] = m_current.y[i] + m_step*(b30*m_k[0][i] + b31*m_k[1][i] + b32*m_k[2][i]);

				m_k[3] = deriv(next.x + a3*m_step, next.y);
				for(size_t i = 0; i < m_num_of_vars; ++i)
					next.y[i] = m_current.y[i] + m_step*(b40*m_k[0][i] + b41*m_k[1][i] + b42*m_k[2][i] + b43*m_k[3][i]);

				m_k[4] = deriv(next.x + a4*m_step, next.y); 
				for(size_t i = 0; i < m_num_of_vars; ++i)
					next.y[i] = m_current.y[i] + m_step*(b50*m_k[0][i] + b51*m_k[1][i] + b52*m_k[2][i] + b53*m_k[3][i] + b54*m_k[4][i]);

				m_k[5] = deriv(next.x + a5*m_step, next.y);
				/*	compute yhat (order 5 RK) and the error directly from m_k values; 
					avoids having to compute order 4 RK explicitly then taking the difference.*/
				double bigerr = 0.0; 
				for(size_t i = 0; i < m_num_of_vars; ++i)
				{
					yhat[i] = m_current.y[i] + m_step*(c0*m_k[0][i] + c2*m_k[2][i] + c3*m_k[3][i] + c4*m_k[4][i] + c5*m_k[5][i]);		
					double err = m_step * fabs(d0*m_k[0][i] + d2*m_k[2][i] + d3*m_k[3][i] + d4*m_k[4][i] + d5*m_k[5][i]);
					//our estimate for the new step should be based on the largest error estimate
					if( fabs(err) > fabs(bigerr) ) bigerr = err;
				}
				hnew = alpha * m_step * sqrt( sqrt( fabs(m_step)*tol/bigerr ) );
				if( fabs(hnew) > 4.0e0 * fabs(m_step) ) hnew = 4.0e0 * m_step;
				if( fabs(hnew) < 1.e-1 * fabs(m_step) ) hnew = 1.e-1 * m_step;
				if( fabs(hnew) > fabs(hmax) )      hnew = hmax;
				if( fabs(hnew) < hmin )
				{
					std::cout << m_step << "->" << hnew << std::endl;
					std::string errmsg = __FUNCTION__;
					errmsg += "Step size too small; exiting integration";
					throw(std::runtime_error(errmsg));
				}

				if( fabs(bigerr) > fabs(m_step*tol) )//repeat integration with new step size
				{
					m_step = hnew; 
				}else // accept the integration at m_current step size
				{
					/*	Reduce the m_step size if it takes the solution beyond the target if set.
						When this occurs we have hit the target so set the bool flag to true in order to store the data externally*/
					if (m_targetState == TARGET && (sign * (next.x + m_step) >= sign * m_nextTarget)){
						m_step = m_nextTarget - next.x;
						m_targetHit = true;
						m_nextTarget += htarg;
					}

					next.x += m_step; 				
					next.y = yhat;
					//adjust next step size to that estimated.
					m_step = hnew;
				}
			}
			m_current = next;
			return next;
		}

		phys::storage::ODE_Storage RKF45::fullSolve(double end)
		{
			phys::storage::ODE_Storage container(m_current);
			state next = m_current;
			yhat = m_current.y;
			while( sign*next.x < sign*end ){
				if (m_targetState == NO_TARGET || m_targetHit){
					container.store(m_current);
					m_targetHit = false;
				}
				while(next.x == m_current.x)//i.e. integrator not advanced
				{
					//compute intermediate values for derivatives
					m_k[0] = deriv(next.x, next.y);
					for(size_t i = 0; i < m_num_of_vars; ++i)
						next.y[i] = m_current.y[i] + b10*m_step*m_k[0][i];

					m_k[1] = deriv(next.x + a1*m_step, next.y);
					for(size_t i = 0; i < m_num_of_vars; ++i)
						next.y[i] = m_current.y[i] + m_step*(b20 * m_k[0][i] + b21*m_k[1][i]);

					m_k[2] = deriv(next.x + a2*m_step, next.y);
					for(size_t i = 0; i < m_num_of_vars; ++i)
						next.y[i] = m_current.y[i] + m_step*(b30*m_k[0][i] + b31*m_k[1][i] + b32*m_k[2][i]);

					m_k[3] = deriv(next.x + a3*m_step, next.y);
					for(size_t i = 0; i < m_num_of_vars; ++i)
						next.y[i] = m_current.y[i] + m_step*(b40*m_k[0][i] + b41*m_k[1][i] + b42*m_k[2][i] + b43*m_k[3][i]);

					m_k[4] = deriv(next.x + a4*m_step, next.y); 
					for(size_t i = 0; i < m_num_of_vars; ++i)
						next.y[i] = m_current.y[i] + m_step*(b50*m_k[0][i] + b51*m_k[1][i] + b52*m_k[2][i] + b53*m_k[3][i] + b54*m_k[4][i]);

					m_k[5] = deriv(next.x + a5*m_step, next.y);
					/*	compute yhat (order 5 RK) and the error directly from m_k values; 
						avoids having to compute order 4 RK explicitly then taking the difference.*/
					double bigerr = 0.0; 
					for(size_t i = 0; i < m_num_of_vars; ++i)
					{
						yhat[i] = m_current.y[i] + m_step*(c0*m_k[0][i] + c2*m_k[2][i] + c3*m_k[3][i] + c4*m_k[4][i] + c5*m_k[5][i]);		
						double err = m_step * fabs(d0*m_k[0][i] + d2*m_k[2][i] + d3*m_k[3][i] + d4*m_k[4][i] + d5*m_k[5][i]);
						//our estimate for the new step should be based on the largest error estimate
						if( fabs(err) > fabs(bigerr) ) bigerr = err;
					}
					hnew = alpha * m_step * sqrt( sqrt( fabs(m_step)*tol/bigerr ) );
					if( fabs(hnew) > 4.0e0 * fabs(m_step) ) hnew = 4.0e0 * m_step;
					if( fabs(hnew) < 1.e-1 * fabs(m_step) ) hnew = 1.e-1 * m_step;
					if( fabs(hnew) > fabs(hmax) )      hnew = hmax;
					if( fabs(hnew) < hmin )
					{
						std::cout << m_step << "->" << hnew << std::endl;
						std::string errmsg = __FUNCTION__;
						errmsg += " Step size too small; exiting integration";
						throw(std::runtime_error(errmsg));
					}

					if( fabs(bigerr) > fabs(m_step*tol) )//repeat integration with new step size
					{
						m_step = hnew; 
					}else // accept the integration at m_current step size
					{
						/*	Reduce the m_step size if it takes the solution beyond the target if set. 
							When this occurs we have hit the target so set the bool flag to true in order to store the data*/
						if ( m_targetState == TARGET && (sign * (next.x + m_step) >= sign * m_nextTarget) ){
							m_step = m_nextTarget - next.x;
							m_targetHit = true;
							m_nextTarget += htarg;
						}

						next.x += m_step; 				
						next.y = yhat;
						//adjust next step size to that estimated.
						m_step = hnew; 
					}
				}//adaptive loop
				m_current = next;
			} //full loop
			container.store(m_current);
			return container;
		}

		phys::storage::ODE_Storage RKF45::fullSolveWrapped(double end)
		{
			phys::storage::ODE_Storage container(m_current);
			do{
				this->wrap();
				if (m_targetState == NO_TARGET || m_targetHit == true){
					container.store(m_current);
				}
				this->solve();
			}while(m_current.x < end);
			this->wrap();
			container.store(m_current); 
			return container;
		}

		/* --- Numerov --- */
		void Numerov::init_method()
		{
			if (m_ptr_diff_eqn->get_order() != 2)
				__error_order();
			if (m_current.num_dims != 1)
				__error_dimensions();
			//set constant used in the algorithm
			m_h_sqr_12 = m_step*m_step / 12.;
			//compute the first two q(x) values. 
			m_k[0] = deriv_B(m_current.x);
			m_k[0] = deriv(m_current.x + m_step, m_k[0]);
		}

		/* As Numerov is designed to only deal with one dimension we know that the vectors
		involved will only ever be of size 2 holding the m_current two values of u(x) or q(x)*/
		state Numerov::solve()
		{
			state next = m_current;
			next.x += m_step;
			next.y[0] = m_current.y[1];
			next.y[1] = 2. * (1. - 5. * m_h_sqr_12 * m_k[0][1]) * m_current.y[1];
			next.y[1] -= (1. + m_h_sqr_12 * m_k[0][0]) * m_current.y[0];
			/*	On entry: m_k[0] = [q(x-m_step)|q(x)] 
				On exit : m_k[0] = [q(x)|q(x+m_step)]	
				where x is the present next.x */
			m_k[0] = deriv(next.x + m_step, m_k[0]); 
			next.y[1] /= 1. + m_h_sqr_12 * m_k[0][1];
			m_current = next;
			return next;
		}

		/*	Set up the storage as a first order ode with one dimension --
			avoids saving repeat values from the previous iteration
		*/ 
		phys::storage::ODE_Storage Numerov::fullSolve(double end)
		{
			size_t N = num_steps(m_current.x, end), count = 0;
			state temp (0.0, 0.0); //arbritrary double values.
			phys::storage::ODE_Storage container(temp); //this just sets up the storage, it doesn't store data.
			state next = m_current;
			while( count++ < N )
			{
				container.store(m_current);
				next.x += m_step;
				next.y[0] = m_current.y[1];
				next.y[1] = 2. * (1. - 5. * m_h_sqr_12 * m_k[0][1]) * m_current.y[1];
				next.y[1] -= (1. + m_h_sqr_12 * m_k[0][0]) * m_current.y[0];
				/*	On entry: m_k[0] = [q(x-m_step)|q(x)] 
					On exit : m_k[0] = [q(x)|q(x+m_step)]	
				where x is the present next.x */
				m_k[0] = deriv(next.x + m_step, m_k[0]); 
				next.y[1] /= 1. + m_h_sqr_12 * m_k[0][1];
				m_current = next;
			}
			container.store(m_current);
			return container;
		}

		void Numerov::__error_order() const
		{
			std::string errmsg = "The Numerov method is specifically used for 2nd order differential\n";
			errmsg += "equations of the form: u''(x) + q(x)u(x) = 0. Please check your code.";
			throw (std::runtime_error(errmsg));
		}
		void Numerov::__error_dimensions() const 
		{
			std::string errmsg = "The Numerov method is use in the one dimensional case only\n";
			errmsg += "Please check your code.";
			throw (std::runtime_error(errmsg));
		}


		/* --- Numerov with source --- */
		void Numerov_s::init_method()
		{
			if (m_ptr_diff_eqn->get_order() != 2) 
				__error_order();
			if (m_current.num_dims != 1) 
				__error_dimensions();
			//set constant used in the algorithm
			m_h_sqr_12 = m_step*m_step / 12.;
			//compute the first two q(x) values. 
			m_k[0] = deriv_B(m_current.x);
			m_k[0] = deriv(m_current.x + m_step, m_k[0]);
			//compute the first two source function values
			m_s[0] = m_src_func(m_current.x);
			m_s[1] = m_src_func(m_current.x + m_step);
		}

		state Numerov_s::solve()
		{
			state next = m_current;
			next.x += m_step;
			m_s[2] = m_src_func(next.x + m_step);
			next.y[0] = m_current.y[1];
			next.y[1] = 2. * (1. - 5. * m_h_sqr_12 * m_k[0][1]) * m_current.y[1];
			next.y[1] -= (1+m_h_sqr_12*m_k[0][0])*m_current.y[0];
			next.y[1] += m_h_sqr_12 * (m_s[2] + m_s[0] + 10. * m_s[1]);
			/* shuffle source function values */
			m_s[0] = m_s[1];
			m_s[1] = m_s[2];

			m_k[0] = deriv(next.x + m_step, m_k[0]); 

			next.y[1] /= 1 + m_h_sqr_12*m_k[0][1];
			m_current = next;
			return next;
		}

		phys::storage::ODE_Storage Numerov_s::fullSolve(double end)
		{
			size_t N = num_steps(m_current.x, end), count = 0;
			phys::storage::ODE_Storage container(m_current);
			state next = m_current;
			while( count++ < N )
			{
				container.store(m_current); 
				next.x += m_step;
				m_s[2] = m_src_func(next.x + m_step);
				next.y[0] = m_current.y[1];
				next.y[1] = 2. * (1. - 5. * m_h_sqr_12 * m_k[0][1]) * m_current.y[1];
				next.y[1] -= (1. + m_h_sqr_12 * m_k[0][0]) * m_current.y[0];
				next.y[1] += m_h_sqr_12 * (m_s[2] + m_s[0] + 10. * m_s[1]);
				m_s[0] = m_s[1];
				m_s[1] = m_s[2];
				m_k[0] = deriv(next.x + m_step, m_k[0]); 
				next.y[1] /= 1 + m_h_sqr_12 * m_k[0][1];
				m_current = next;
			}
			container.store(m_current);
			return container;
		}

		void Numerov_s::__error_order() const
		{
			std::string errmsg = "The Numerov_s method is specifically used for 2nd order differential\n";
			errmsg += "equations of the form: u''(x) + q(x)u(x) = s(x). Please check your code.";
			throw (std::runtime_error(errmsg));
		}
		void Numerov_s::__error_dimensions() const
		{
			std::string errmsg = "The Numerov_s method is use in the one dimensional case only\n";
			errmsg += "Please check your code.";
			throw (std::runtime_error(errmsg));
		}

		/* ---------General error messages for the ODE solver class --------------------------------*/

		void ODE_solver::__error_message_order_1() const
		{
			std::string errmsg = "You have called/defined an order one ODE but have\n";
			errmsg += "specified first derivative varibles which are NOT required.\n";
			errmsg += "Please check your code.";
			std::cout << errmsg << std::endl;
			throw(std::runtime_error(errmsg));
		}

		void ODE_solver::__error_message_order_2() const
		{
			std::string errmsg = "You have called/defined an order two ODE but have\n";
			errmsg += "not specified first derivative varibles which ARE required.\n";
			errmsg += "Please check your code.";
			std::cout << errmsg << std::endl;
			throw(std::runtime_error(errmsg));
		}
	}
}
