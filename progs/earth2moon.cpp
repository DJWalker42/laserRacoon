#include <ODESolvers.h>
#include <Storage.h>
#include <Visualise.h>
#include <PhysicalUnits.h>
#include <opencv2/opencv.hpp>

/*
	***Program to Compute a possible journey of a space craft to the moon from low earth orbit***

	Uses the Runge-Kutta-Fehlberg algorithm to compute the craft's trajactory under the influence 
	of the Earth's and Moon's gravitation fields (other gravitational influences are assumed insignificant).
	It is assumed the craft fires some booster rockets which is considered an instantaneous kick to its
	initial velocity. This assumption is okay as the rockets only fire for a few minutes whereas the total
	journery will be several days at least.

	The coordinate origin lies at the centre of mass between the Earth and the Moon. We assume the Moon-Earth system
	moves in a circular orbit about this centre of mass. Motion of the space craft is reistricted to the plane of orbit.
	
	Units are lunar distance, i.e. the moon's radius of orbit, and a sidereal month for time.
*/



/* Global variables*/
double e_coords[2] = {0.0}; 
double m_coords[2] = {0.0}; 

const double MM = phys::constants::moon_mass/phys::constants::earth_mass/1.e2;
const double rm = 1 / ( MM + 1 );
const double re = rm * MM;

double phys::diffs::User_eqn::differential_function(double t, const std::vector<double>& y, int N, int c)
{
	using namespace phys::constants;
	e_coords[0] = re * cos(2*PI*t + PI);
	e_coords[1] = re * sin(2*PI*t + PI); 

	m_coords[0] = rm * cos(2*PI*t);
	m_coords[1] = rm * sin(2*PI*t); 

	double de, dm, de3, dm3;
	de = sqrt( (y[0] - e_coords[0]) * (y[0] - e_coords[0]) + (y[1] - e_coords[1]) * (y[1] - e_coords[1]) );
	dm = sqrt( (y[0] - m_coords[0]) * (y[0] - m_coords[0]) + (y[1] - m_coords[1]) * (y[1] - m_coords[1]) );

	de3 = de*de*de;
	dm3 = dm*dm*dm;

	return ( -GE * ( ((y[c]-e_coords[c])/de3) + (MM * (y[c] - m_coords[c])/dm3) ) );  
}

int main()
{
	using namespace phys::constants;

	int rows = 1000;
	int cols = 1000;

	std::string winname = "Plot"; 

	double t = 0.0;
	std::vector<double> position;
	std::vector<double> velocity;
	double step = 0.1;
	double theta = 3.0; 
	double dv = 3.084e1; //"instantaneous" kick to velocity due to space craft firing rockets.


	//space craft starts out in low earth orbit with an instantaneous kick to its tangential velocity.
	position.push_back( (6.8e-2 * sin(theta)/lunar_dist) - re );					//initial x position
	position.push_back( (6.8e-2 * cos(theta)/lunar_dist) - e_coords[1] );			//initial y position
	velocity.push_back( (7.657e1 + dv) * cos(theta) * sidereal_month/lunar_dist );	//initial x velocity
	velocity.push_back( -(7.657e1 + dv) * sin(theta) * sidereal_month/lunar_dist);	//initial y velocity

	phys::ode::state initial_system(t, position, velocity);

	phys::diffs::Diff_eqn*travel = new phys::diffs::User_eqn;
	travel->set_order(2);

	phys::ode::ODE_solver*rk45 = new phys::ode::RKF45(travel, initial_system, step);

	double moon_radius_check = moon_radius/lunar_dist/1.e2;
	std::cout << "Target: " << moon_radius_check << "\n";

	phys::ode::state next = initial_system;

	double distance_to_moon = sqrt( (next.y[0] - m_coords[0]) * (next.y[0] - m_coords[0]) + 
			(next.y[1] - m_coords[1]) * (next.y[1] - m_coords[1]) );

	cv::Mat BG(rows, cols, CV_8UC3, cv::Scalar(0,0,0));

	double scale = rows/3;

	//coordinate origin == centre of mass for the system: found at (cols/2, rows/2) in the window.

	int i_earth_radius = int(scale*earth_radius/lunar_dist/1.e2);
	int i_moon_radius = int(scale*moon_radius_check);

	cv::Mat display, display2;

	cv::namedWindow(winname);//, CV_WINDOW_KEEPRATIO);

	display2 = BG.clone();

	//main loop: break if space craft is more than two lunar distances from moon or has been travelling longer than one sidereal month.
	while(distance_to_moon < 2.0 && next.x < 1.0) 
	{
		display = BG.clone();
		cv::Point earth_centre( int(e_coords[0] * scale) + cols/2, rows/2 - int(e_coords[1] * scale) );
		cv::Point moon_centre( int(m_coords[0] * scale) + cols/2, rows/2 - int(m_coords[1] * scale) );
		cv::Point space_craft( int(scale*next.y[0]) + cols/2, rows/2 - int(scale*next.y[1]) );

		cv::circle(display, earth_centre, i_earth_radius, cv::Scalar(255, 0, 0), -1);
		cv::circle(display, moon_centre, i_moon_radius, cv::Scalar(255, 255, 255), -1);
		cv::circle(display2, space_craft, 1, cv::Scalar(0,0,255), -1);

		display2.copyTo(display,display2);

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
	getchar();
	return 0;
}
