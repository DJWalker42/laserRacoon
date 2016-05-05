#include <ODESolvers.h>
#include <PhysicalUnits.h>
#include <Visualise.h>

/*
	***PROGRAM TO COMPUTE THE ORBIT OF HALLEY'S COMET***
	
	Uses the adaptive Runge-Kutta-Fehlberg algorithm to compute the orbit of Halley's comet around the Sun.
	Note that Halley's comet is defined by its aphelion position and tangential velocity at that point.
	To keep computer friendly numbers during the computation we use the Astronomical Unit (AU) for the distance 
	unit and the number of seconds in an Earth year (yrs) for the time unit (both found in PhysicalUnits.h).

	Here we assume no other gravitational interaction with other bodies in the solar system. 
*/

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
		viewer.plot(data, Viewer::XY);
	}

	data.write("./data/cometOutput.txt");

	return 0; 
}
