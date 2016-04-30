#include <Storage.h>
#include <ODESolvers.h>
#include <Visualise.h>

int main ()
{
	double t = 0.0;
	double p = 1.0;
	double v = 0.0;
	double step = 0.02;
	double finish = 100;

	phys::ode::state initial_system(t,p,v); 
	phys::diffs::Diff_eqn* eqn = new phys::diffs::Duffing(-1.0, 1.0, 0.2, 0.3, 1.0);
	phys::ode::ODE_solver* rk4_solver = new phys::ode::RK4(eqn, initial_system, step);
	phys::storage::ODE_Storage container = rk4_solver->fullSolve(finish);

	delete eqn;
	delete rk4_solver;

	phys::visual::Viewer viewer;
	viewer.plot( container, phys::visual::Viewer::PHASE );

	//container.write("../Data/output.txt");

	return 0;
}
