#include <Storage.h>
#include <ODESolvers.h>
#include <Visualise.h>

int main ()
{
	double t = 0.0;
	double p = 1.0;
	double v = 0.0;
	double step = 0.02;
	double finish = 10.0;

	phys::ode::state initial_system(t,p,v); 
	phys::diffs::Diff_eqn* shm_eqn = new phys::diffs::SHM_eqn();
	phys::ode::ODE_solver* rk4_solver = new phys::ode::RKF45(shm_eqn, initial_system, step);
	phys::storage::ODE_Storage container = rk4_solver->fullSolve(finish);

	delete shm_eqn;
	delete rk4_solver;

	phys::visual::Viewer viewer;

	viewer.plot( container);//, phys::visual::Viewer::PHASE ); 

	return 0;
}
