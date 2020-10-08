#include <Storage.h>
#include <ODESolvers.h>
#include <Visualise.h>

/*
	***PROGRAM TO COMPUTE THE SOLUTION TO THE DUFFING OSCILLATOR***
	
	The Duffing oscillator shows chaotic motion for given sets of parameters and initial conditions.
	This can be seen in the phase space plot of the solution. 
	See http://www.scholarpedia.org/article/Duffing_oscillator for a technical explanation.
*/

int main ()
{
	double t = 0.0;
	double p = 1.0;
	double v = 0.0;
	double step = 0.02;
	double finish = 100;

	phys::state initial_system(t,p,v);
	phys::Diff_eqn* eqn = new phys::Duffing(-1.0, 1.0, 0.2, 0.3, 1.0);
	phys::ODESolver* rk4_solver = new phys::RK4(eqn, initial_system, step);
	phys::ODEStorage container = rk4_solver->fullSolve(finish);

	delete eqn;
	delete rk4_solver;

	phys::Viewer viewer;
	viewer.plot( container, phys::Viewer::PHASE );

	//container.write("../Data/output.txt");

	return 0;
}
