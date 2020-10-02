#include <Storage.h>
#include <ODESolvers.h>
#include <Visualise.h>

/*
	*** SIMULATION OF SIMPLE HARMONIC MOTION ***
	
	This program uses the Runge-Kutta algorithm to simulate SHM for a given initial state.

	We use base class pointers with new and delete to show you the syntax; in modern C++
	you'd use smart pointers. We do this as a educational point. You could instead just
	create direct objects of type SHM and RKF45, with appropriate use of the "address of"
	operator (&).
*/


int main ()
{
	double t = 0.0; //start time
	double p = 1.0; //position
	double v = 0.0; //velocity
	double step = 0.02; //integration step
	double finish = 10.0; //end time

	phys::ode::state initial_system(t,p,v); 
	phys::diffs::Diff_eqn* shm_eqn = new phys::diffs::SHM();
	phys::ode::ODESolver* rk4_solver = new phys::ode::RKF45(shm_eqn, initial_system, step);
	phys::storage::ODEStorage container = rk4_solver->fullSolve(finish);

	delete shm_eqn;
	delete rk4_solver;

	phys::visual::Viewer viewer;

	//uncomment the second argument to see the phase-space plot
	viewer.plot(container);//, phys::visual::Viewer::PHASE ); 

	return 0;
}
