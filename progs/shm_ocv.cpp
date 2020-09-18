#include <Storage.h>
#include <ODESolvers.h>
#include <Visualise.h>

/*
	*** SIMULATION OF SIMPLE HARMONIC MOTION ***
	
	This program uses the Runge-Kutta algorithm to simulate SHM for a given initial state.
*/


int main ()
{
	double t = 0.0; //start time
	double p = 1.0; //position
	double v = 0.0; //velocity
	double step = 0.02; //integration step
	double finish = 10.0; //end time

	phys::ode::state initial_system(t,p,v); 
	phys::diffs::Diff_eqn* shm_eqn = new phys::diffs::SHM_eqn();
	phys::ode::ODE_solver* rk4_solver = new phys::ode::RKF45(shm_eqn, initial_system, step);
	phys::storage::ODE_Storage container = rk4_solver->fullSolve(finish);

	delete shm_eqn;
	delete rk4_solver;

	phys::visual::Viewer viewer;

	//uncomment the second argument to see the phase-space plot
	viewer.plot(container);//, phys::visual::Viewer::PHASE ); 

	return 0;
}
