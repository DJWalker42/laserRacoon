#include <ODESolvers.h>
#include <PhysicalUnits.h>
#include <Visualise.h>


int main()
{
	phys::stdVec_d position(2), velocity(2);
	double step = .1, t = .0; 

	position[0] = 0.;
	position[1] = 35.1; //distance of comet at aphelion in AU
	velocity[0] = 8.79e2 * yrs * 1.e-4 / AU;   //speed of comet at aphelion in AU per year 
	velocity[1] = 0.;

	phys::ode::state initial(t, position, velocity); 

	phys::diffs::Gravity grav;
	phys::ode::RKF45 rkf45(&grav, initial, step);

	phys::storage::ODE_Storage data = rkf45.fullSolve(300.);

	{
		using namespace phys::visual;
		Viewer viewer;
		viewer.plot(data);//, Viewer::XY);
	}

	data.write("./data/cometOutput.txt");

	return 0; 
}
