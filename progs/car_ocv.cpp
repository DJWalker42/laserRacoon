#include <iostream>
#include <vector>
#include <Visualise.h>

class car{
	double speed;
	double mass;
	double power;		//horse power - maybe thrust?
	double drag_coeff;
public:
	car() :	speed(0.),
			mass(0.),
			power(0.), 
			drag_coeff(0.) {}

	car(	double m, 
			double p,
			double dc) :
			speed(0.),
			mass(m),
			power(p),
			drag_coeff(dc) {}
	
	double accelerate( double amount );
	double brake( double amount ); 
};


double car::accelerate( double amt )
{
	speed += amt*(power - drag_coeff*speed)/mass; 
	return speed;
}

double car::brake( double amt )
{
	if(speed > 0.)
		speed -= amt*(power + drag_coeff*speed)/mass; 
	return speed > 0. ? speed : 0;
}

int main()
{
	car banger(10., 4., 3.);
	car sports(5., 10., 1.);

	phys::stdVec_d time;
	std::vector<phys::stdVec_d> velocity(2);

	double amt = .5;

	velocity[0].push_back(0.);
	velocity[1].push_back(0.);
	time.push_back(0.);

	for(size_t i = 1; i < 100; ++i)
	{
		velocity[0].push_back(banger.accelerate(amt));
		velocity[1].push_back(sports.accelerate(amt));
		time.push_back(i);
	}

	amt = .1;
	for(size_t i = 100; i < 200; ++i)
	{
		velocity[0].push_back(banger.brake(amt));
		velocity[1].push_back(sports.brake(amt));
		time.push_back(i);
	}

	phys::Viewer viewer;
	viewer.set_x_name("time/s");
	viewer.set_y_name("velocity/m/s");
	viewer.plot(time, velocity);

	return 0;
}
