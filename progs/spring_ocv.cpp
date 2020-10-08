#include <Visualise.h>

const double PI = 4. * atan(1.);

class Spring{
private:
	double _position;
	double _velocity;
	double _springConst;
	double _dragCoeff;
public:
	Spring(	double position, 
			double velocity,
			double springConstant, 
			double dragCoefficient) :
			_position(position),
			_velocity(velocity),
			_springConst(springConstant),
			_dragCoeff(dragCoefficient)
	{}

public:
	void update( double step );

	double get_pos() const {return _position;}
	double get_vel() const {return _velocity;}
};

//Forward Euler for a 2nd order differential
void Spring::update( double h )
{
	double old_pos = _position;
	_position += _velocity * h;
	_velocity -= _springConst * old_pos * h + _dragCoeff * _velocity * h;
}

int main()
{	
	double initPos = 1.;
	double initVel = 0.;
	double springConst = 4. * PI * PI; //this value gives an undamped time period of 1 unit (seconds)
	double dragCoeff = .0;

	std::cout << "Specify a drag coefficient (double): ";
	std::cin >> dragCoeff;

	//FIXME: validate input please, dragCoeff could be anything!!

	//construct a spring object with the initial condition, and the physical constants
	Spring spring(initPos, initVel, springConst, dragCoeff); 

	//names for the "state" of the spring
	std::vector<std::string> stateID { "Pos.", "Vel." };

	//storage for the data using the names specified. This allows us to store time, position, and velocity
	phys::Storage<double> state("time", stateID);

	//initialise a time step
	double step = 0.00001;
	for(double t = 0.; t < 10.; t += step)
	{
		state.store(t, spring.get_pos(), spring.get_vel());
		spring.update(step);
	}

	//create a viewer using the default constructor
	phys::Viewer viewer;
	viewer.withLines(); //option to draw lines between points
	viewer.plot(state, true); //plot the data using split vertical axes

	//press any key to advance the data

	return 0;
}
