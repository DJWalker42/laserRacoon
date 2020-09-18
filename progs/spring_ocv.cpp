#include <Visualise.h>

const double PI = 4. * atan(1.);

class spring{
public:
	spring(	double position, 
			double velocity,
			double springConstant, 
			double dragCoefficient) :
			m_position(position),
			m_velocity(velocity),
			m_springConst(springConstant),
			m_dragCoeff(dragCoefficient)
	{}

	void update( double step );

	double get_pos() const {return m_position;}
	double get_vel() const {return m_velocity;}

private:
	double m_position;
	double m_velocity;
	double m_springConst;
	double m_dragCoeff;
};

void spring::update( double h )
{
	double old_pos = m_position;
	m_position += m_velocity * h;
	m_velocity -= m_springConst * old_pos * h + m_dragCoeff * m_velocity * h; 
}

int main()
{	
	double initPos = 1.;
	double initVel = 0.;
	double springConst = 4. * PI * PI; //this value gives an undamped time period of 1 unit (seconds)
	double dragCoeff = .0;

	std::cout << "Specify a drag coefficient (double): ";
	std::cin >> dragCoeff;

	//construct a spring object with the initial condition, and the physical constants
	spring aSpring(initPos, initVel, springConst, dragCoeff); 

	//names for the "state" of the spring
	std::vector<std::string> stateID = { "Pos.", "Vel." };

	//storage for the data using the names specified. This allows us to store time, position, and velocity
	phys::storage::Storage<double> state("time", stateID);

	//initialise a time step
	double step = 0.00001;
	for(double t = 0.; t < 10.; t += step)
	{
		state.store(t, aSpring.get_pos(), aSpring.get_vel());
		aSpring.update(step);
	}

	//create a viewer using the default constructor
	phys::visual::Viewer viewer;
	viewer.withLines(); //option to draw lines between points
	viewer.plot(state, true); //plot the data using split vertical axes

	return 0;
}
