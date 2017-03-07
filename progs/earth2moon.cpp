#include <ODESolvers.h>
#include <Storage.h>
#include <Visualise.h>
#include <PhysicalUnits.h>
#include <opencv2/opencv.hpp>

/*
	***Program to Compute a possible journey of a space craft to the moon from low earth orbit***

	Uses the Runge-Kutta-Fehlberg algorithm to compute the craft's trajectory under the influence
	of the Earth's and Moon's gravitation fields (other gravitational influences are assumed insignificant).
	It is assumed the craft fires some booster rockets that provides an instantaneous kick to its
	initial velocity. This assumption is okay as the rockets only fire for a few minutes whereas the total
	journey will be several days at least.

	The coordinate origin lies at the centre of mass between the Earth and the Moon. We assume the
	Moon-Earth system moves in a circular orbit about this centre of mass. Motion of the space craft
	is restricted to the plane of orbit. We assume the initial position of the Earth and Moon is on
	the horizontal axis of our coordinate system. This is okay as we can choose any arbitrary starting
	angular position for our space craft. The craft moves in a clockwise orbit as we look down on the
	x-y plane (this can be changed by swapping the parity of the velocity components), whereas the
	earth-moon system moves in an anti-clockwise sense.
	
	Units are lunar distance (Earth to moon distance), Earth's mass, and a sidereal month for time.

	Improvements:
	Is the assumption that the Earth-Moon orbit is circular a valid one?

	Are the other gravitational sources in the Solar System really insignificant?

	Are sidereal months the best unit of time for this computation? How about (sidereal) days instead?

	How much of difference does it make if we accurately model the booster rockets rather than assuming
	they are instantaneous?

	We have "Success" if the space craft distance from the (centre of the) Moon is less than its radius,
	what does this imply? Instead, could we get the craft into a stable orbit with the Moon?

	Assuming we've managed to land intact on the Moon, could you get the craft, and by implication, the
	astronauts back to Earth?
*/


/* Global variables*/
const double MM = phys::constants::moon_mass/phys::constants::earth_mass/1.e2;
const double rm = 1 / ( MM + 1 ); //distance of the centre of the Moon from the centre of mass
const double re = rm * MM; //distance of the centre of the Earth from the centre of mass

//in our coordinate system assume the centres of the Earth and Moon lie on the horizontal axis.
double e_coords[2] = {-re, 0.0};
double m_coords[2] = {rm, 0.0};

// Craft's low Earth orbit radius from the centre of the Earth, and corresponding tangential speed.
// Note the craft's orbit is assumed circular.
const double orbit_radius = 6.8; //units: km
const double orbit_speed = 7.657; //units: km/s


//define our differential_function that is called during the integration
//parameters are time, craft x-y position, number of dimensions in system (2 in this case),
//and the specific dimension to currently compute.
//Note that 'y' also contains the craft's x-y velocity in indices 2 and 3, respectively, but these are unused.

double phys::diffs::User_eqn::differential_function(double t, const std::vector<double>& y, int N, int c)
{
	using namespace phys::constants; // for PI and GE

	//Compute where the Earth and Moon have moved to in their circular orbit about centre-of-mass
	e_coords[0] = re * cos(2*PI*t + PI);
	e_coords[1] = re * sin(2*PI*t + PI); 
	m_coords[0] = rm * cos(2*PI*t);
	m_coords[1] = rm * sin(2*PI*t); 

	double de, dm; //distance of craft from earth and moon respectively
	de = sqrt( (y[0] - e_coords[0]) * (y[0] - e_coords[0]) + (y[1] - e_coords[1]) * (y[1] - e_coords[1]) );
	dm = sqrt( (y[0] - m_coords[0]) * (y[0] - m_coords[0]) + (y[1] - m_coords[1]) * (y[1] - m_coords[1]) );

	//Cube those values for force (acceleration) calculation
	double de3 = de*de*de;
	double dm3 = dm*dm*dm;

	return ( -GE * ( ((y[c]-e_coords[c])/de3) + (MM * (y[c] - m_coords[c])/dm3) ) );  
}

int main()
{
	using namespace phys::constants;

	int rows = 1000;
	int cols = 1000;

	std::string winname = "Plot"; 

	double initial_time = 0.0;
	std::vector<double> position;
	std::vector<double> velocity;
	double step = 0.1;
	double theta; //Space craft's angular position w.r.t. the centre of mass of the earth-moon system
	double dv; //"instantaneous" kick to velocity due to space craft firing rockets.

	//Get user to supply theta and dv. Note theta is in radians. You can put in any boost value you
	//like but for sensible results make it less than 5.0. To check that the low earth orbit is
	//stable make the boost zero.

	std::cout << "Space craft's angular position (rads.): ";
	std::cin >> theta;
	std::cout << "Amount of velocity gained from boost [0.0, 5.0] (km/s): ";
	std::cin >> dv;
	dv *= 10.; //for computational units


	//Space craft starts out in low earth orbit with an instantaneous kick to its tangential velocity.
	//The magic numbers 1.e-2 and 1.e1 are required to get the correct scale for our choice of units found in Physical_Units.cpp
	//In current parity (+vx, -vy) orbit is clockwise, swap to get an anti-clockwise orbit (-vx, +vy).
	position.push_back( (orbit_radius * 1.e-2 * sin(theta)/lunar_dist) - re );					//initial x position
	position.push_back( (orbit_radius * 1.e-2 * cos(theta)/lunar_dist) - e_coords[1] );			//initial y position
	velocity.push_back( (orbit_speed * 1.e1 + dv) * cos(theta) * sidereal_month/lunar_dist );	//initial x velocity
	velocity.push_back( -(orbit_speed * 1.e1 + dv) * sin(theta) * sidereal_month/lunar_dist);	//initial y velocity

	phys::ode::state initial_system(initial_time, position, velocity);

	phys::diffs::Diff_eqn *travel = new phys::diffs::User_eqn;
	travel->set_order(2);

	phys::ode::ODE_solver *rk45 = new phys::ode::RKF45(travel, initial_system, step);

	double moon_radius_check = moon_radius/lunar_dist/1.e2;
	std::cout << "Target: " << moon_radius_check << "\n";

	phys::ode::state next = initial_system;

	double distance_to_moon = sqrt( (next.y[0] - m_coords[0]) * (next.y[0] - m_coords[0]) + 
			(next.y[1] - m_coords[1]) * (next.y[1] - m_coords[1]) );

	cv::Mat BG(rows, cols, CV_8UC3, cv::Scalar(0,0,0)); //black background image

	double scale = rows/3; //ensure earth-moon motion remains in confines of the image

	//coordinate origin == centre of mass for the system: found at (cols/2, rows/2) in the window.

	int i_earth_radius = int(scale*earth_radius/lunar_dist/1.e2);
	int i_moon_radius = int(scale*moon_radius_check);

	cv::Mat display, display2;

	cv::namedWindow(winname);//, CV_WINDOW_KEEPRATIO);

	display2 = BG.clone(); //so we can see the trajectory of the craft

	//main loop: break if space craft is more than five lunar distances from moon or has been travelling longer than one sidereal month.
	while(distance_to_moon < 5.0 && next.x < 1.0)
	{
		display = BG.clone();
		cv::Point earth_centre( int(e_coords[0] * scale) + cols/2, rows/2 - int(e_coords[1] * scale) );
		cv::Point moon_centre( int(m_coords[0] * scale) + cols/2, rows/2 - int(m_coords[1] * scale) );
		cv::Point space_craft( int(scale*next.y[0]) + cols/2, rows/2 - int(scale*next.y[1]) );

		cv::circle(display, earth_centre, i_earth_radius, cv::Scalar(255, 0, 0), -1);	//earth is blue, obvs.
		cv::circle(display, moon_centre, i_moon_radius, cv::Scalar(255, 255, 255), -1);	//moon is white, obvs.
		cv::circle(display2, space_craft, 1, cv::Scalar(0,0,255), -1); //space craft is red, not obvs. but is clear.

		display2.copyTo(display,display2); //copies craft trajectory to our "main" display image.

		next = rk45->solve();

		distance_to_moon = sqrt( (next.y[0] - m_coords[0]) * (next.y[0] - m_coords[0]) + 
			(next.y[1] - m_coords[1]) * (next.y[1] - m_coords[1]) );

		std::cout << "Current distance to moon: " << distance_to_moon << "\n";
		if(distance_to_moon <= moon_radius_check)  
		{
			cv::putText( display, "Success!", moon_centre + cv::Point(-40,-40), cv::FONT_HERSHEY_PLAIN, 3.0, cv::Scalar(127,255,60), 3 );
			std::cout << "Success!\n" ;
			std::cout << "time taken = " << next.x*sidereal_month*1e6/60/60/24 << " days" << std::endl;
			cv::imshow(winname, display);
			cv::waitKey();
			exit(EXIT_SUCCESS); //break out here on success
		}

		cv::imshow(winname, display); 
		char key = cv::waitKey(100);
		switch(key){
		case 27:
			cv::destroyAllWindows();
			return 0; //break out on an Esc press
		}
	}

	//if here while loop conditional returned false, display a failure message.
	std::cout << "Missed!\n";
	std::cout << "You are now " << distance_to_moon << " lunar distances from the moon\n";
	std::cout << "You have been travelling for " << next.x*sidereal_month*1e6/60/60/24 << " days" << std::endl;
	cv::waitKey();
	return 0;
}
