#include <ODESolvers.h>
#include <Storage.h>
#include <Visualise.h>
#include <PhysicalUnits.h>
#include <opencv2/opencv.hpp>

/* Global variables*/
double e_coords[2] = {0.0};
double m_coords[2] = {0.0}; 

const double MM = moon_mass/earth_mass/1.e2;
const double rm = 1 / ( MM + 1 );
const double re = rm * MM;

double phys::diffs::User_eqn::differential_function(double t, const std::vector<double>& y, int N, int c)
{
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
	int rows = 1000;
	int cols = 1000;

	std::string winname = "Plot"; 

	double t = 0.0;
	std::vector<double> position;
	std::vector<double> velocity;
	double step = 0.1;
	double theta = 3.0; 
	double dv = 3.084e1;

	position.push_back( (6.8e-2 * sin(theta)/lunar_dist) - re );					//initial x position
	position.push_back( (6.8e-2 * cos(theta)/lunar_dist) - e_coords[1] );			//initial y position
	velocity.push_back( (7.657e1 + dv) * cos(theta) * sidereal_month/lunar_dist );	//initial x velocity
	velocity.push_back( -(7.657e1 + dv) * sin(theta) * sidereal_month/lunar_dist);	//initial y velocity

	phys::ode::state initial_system(t, position, velocity);

	phys::diffs::Diff_eqn*travel = new phys::diffs::User_eqn;
	travel->set_order(2);

	phys::ode::ODE_solver*rk45 = new phys::ode::RKF45(travel, initial_system, step);

	double m_r_chk = moon_radius/lunar_dist/1.e2;
	std::cout << "Target: " << m_r_chk << "\n";

	phys::ode::state next = initial_system;

	double d_to_m = sqrt( (next.y[0] - m_coords[0]) * (next.y[0] - m_coords[0]) + 
			(next.y[1] - m_coords[1]) * (next.y[1] - m_coords[1]) );

	cv::Mat BG(rows, cols, CV_8UC3, cv::Scalar(0,0,0));

	double scale = rows/3;

	//origin == centre of mass: found at (cols/2, rows/2)


	int e_radius = int(scale*earth_radius/lunar_dist/1.e2);
	int m_radius = int(scale*m_r_chk);

	cv::Mat display, display2;

	cv::namedWindow(winname);//, CV_WINDOW_KEEPRATIO);

	display2 = BG.clone();

	while(d_to_m < 2.0 && next.x < 1.0) 
	{
		display = BG.clone();
		cv::Point eh( int(e_coords[0] * scale) + cols/2, rows/2 - int(e_coords[1] * scale) );
		cv::Point mn( int(m_coords[0] * scale) + cols/2, rows/2 - int(m_coords[1] * scale) );
		cv::Point sp( int(scale*next.y[0]) + cols/2, rows/2 - int(scale*next.y[1]) );

		cv::circle(display, eh, e_radius, cv::Scalar(255,0,0), -1);
		cv::circle(display, mn, m_radius, cv::Scalar(255,255,255), -1);
		cv::circle(display2, sp, 1, cv::Scalar(0,0,255), -1);

		display2.copyTo(display,display2);

		next = rk45->solve();

		d_to_m = sqrt( (next.y[0] - m_coords[0]) * (next.y[0] - m_coords[0]) + 
			(next.y[1] - m_coords[1]) * (next.y[1] - m_coords[1]) );

		std::cout << d_to_m << "\n";
		if(d_to_m <= m_r_chk)  
		{
			cv::putText( display, "Success!", mn + cv::Point(-40,-40), cv::FONT_HERSHEY_PLAIN, 3.0, cv::Scalar(127,255,60), 3 );
			std::cout << "Success!\n" ;
			std::cout << "time taken = " << next.x*sidereal_month*1e6/60/60/24 << " days" << std::endl;
			cv::imshow(winname, display);
			cv::waitKey();
			exit(EXIT_SUCCESS);
		}

		cv::imshow(winname, display); 
		char key = cv::waitKey(100);
		switch(key){
		case 27:
			cv::destroyAllWindows();
			return 0;
		}
	}

	std::cout << "Missed!\n";
	std::cout << "You are now " << d_to_m << " lunar distances from the moon\n";
	std::cout << "You have been travelling for " << next.x*sidereal_month*1e6/60/60/24 << " days" << std::endl;
	getchar();
	return 0;
}
