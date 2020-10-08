#include <ODESolvers.h>
#include <Visualise.h>
#include <sstream>

#include <opencv2/imgproc.hpp>

/*
	***Program to simulate a double pendulum***

	Uses the adaptive Runge-Kutta-Fehlberg algorithm to compute the motion of a double pendulum system.
	Note that drag has not been considered (the system is undamped) and the motion is constrained to a
	plane only.

	The data is then animated using OpenCV. After quitting the animation (q,Q, or Esc) a phase-space
	plot is shown. For the current parameters: pendulum_1 length 1 mass 1, pendulum_2 length .5 mass 2,
	and run time of 100 this plot looks interesting.
*/

const int rows = 500;
const int cols = rows;

cv::Point origin( cols/2, rows/2 ); 

const int end = 100; //final time
const double scale = 100.; //pixels per unit
const int delay = 30; //milliseconds delay between animation frames. To have animation run in "real time" make the time step = delay / 1000.

std::string winname = "Double Pendulum";


double l[] = {1., .5}; //basic array representing the lengths of the two pendula
double m[] = {1., 2.}; //basic array representing the masses of the two pendula


//here the parameter 'y' contains the angular positions and velocities for both pendula.

double phys::User_eqn::differential_function(double t, const phys::stdVec_d& y, int N, int c)
{
	double G = 9.81;
	double del = y[1] - y[0];
	double denom1 = (m[0] + m[1])*l[0] - m[1]*l[0]*cos(del)*cos(del);
	double denom2 = (l[1]/l[0]) * denom1;
	double retval = 0.;
	switch(c) // c will only ever be zero or one
	{
	case 0:
		retval = m[1]* l[0] * y[2] * y[2] * sin(del) * cos(del);
		retval += m[1] * G * sin(y[1]) * cos(del);
		retval += m[1] * l[1] * y[3] * y[3] * sin(del) ;
		retval -= (m[0] + m[1]) * G * sin( y[0]);
		retval /= denom1;
		break;
	case 1:
		retval = -m[1] * l[1] * y[3] * y[3] * sin(del) * cos(del);
		retval += (m[0] + m[1]) * G * sin(y[0]) * cos(del); 
		retval -= (m[0] + m[1]) * l[0] * y[2] * y[2] * sin(del);
		retval -= (m[0] + m[1]) * G * sin(y[1]); 
		retval /= denom2;
		break;
	}
	return retval;
}


bool show(const cv::Mat& frame)
{
	cv::imshow(winname, frame);
	char key = char(cv::waitKey(delay)); //delay to match the time step used in the RK4
	switch (key)
	{
	case 'q':
	case 'Q': //fall through wanted
	case 27: //27 == Esc key
		cv::destroyWindow(winname);
		return false;
	}
	return true;
}

void animate( const phys::ODEStorage& data )
{
	phys::stdVec_d theta_1 = data.get_dependent(1);
	phys::stdVec_d theta_2 = data.get_dependent(2);

	size_t n = theta_1.size(); // == theta_2.size()

	std::vector<cv::Point> p1(n);
	std::vector<cv::Point> p2(n);

	for (size_t i = 0; i < n; ++i)
	{
		p1[i] = cv::Point(origin.x + int(scale * l[0] * sin(theta_1[i])), 
			int(scale * l[0] * cos(theta_1[i])) + origin.y);

		p2[i] = cv::Point(int(scale * l[1] * sin(theta_2[i])) + p1[i].x, 
			int(scale * l[1] * cos(theta_2[i])) + p1[i].y );

	}

	cv::Mat frame0(rows, cols, CV_8U, cv::Scalar(0));
	cv::Mat frame;

	size_t frame_count = 0;

	do
	{
		//computes time in seconds - used to display the time in the window
		double time = static_cast<double>( (frame_count * delay)%(end * 1000) )/ 1000.;

		std::stringstream ss;
		ss << time << " s"; 

		frame = frame0.clone();
		size_t idx = frame_count++ % n;

		//this displays the time in the top left corner
		cv::putText(frame, ss.str(), cv::Point(20,20), cv::FONT_HERSHEY_PLAIN, 1., cv::Scalar(255)); 

		//ends of the pendula.
		cv::Point end1 = cv::Point(p1[idx].x, p1[idx].y); 
		cv::Point end2 = cv::Point(p2[idx].x, p2[idx].y);

		cv::line(frame, origin, end1, cv::Scalar(255), 2);  //white
		cv::line(frame, end1, end2, cv::Scalar(127), 2);	//grey
	}while( show(frame) );
}


int main()
{
	//theta are the initial pendulum angles in radians
	//omega are the initial pendulum angular velocities in radians per second

	phys::stdVec_d theta(2), omega(2);

	//double rad = phys::constants::PI/180;
	theta[0] = 0.;//120. * rad;
	theta[1] = 0.;//-10. * rad;
	omega[0] = 2.;//0.;
	omega[1] = 7.;//0.;

	double t = 0.;
	double step = static_cast<double>(delay) / 1000.;

	phys::User_eqn pendulum(2); //2nd order ODE
	phys::state initial(t, theta, omega);

	//by using a target we save data at the target steps only (the solver still adapts in-between)
	phys::RKF45 rkf45(&pendulum, initial, step, phys::TARGET);

	//wraps the solution to [-pi, +pi]
	phys::ODEStorage data = rkf45.fullSolveWrapped(static_cast<double>(end));

	animate(data);

	phys::Viewer viewer;
	viewer.set_plot_name("Foxy?"); 
	viewer.withLines();
	viewer.plot(data, phys::Viewer::PHASE, 1); //the final argument specifies that we want to see the data for the *second* pendulum.

	return 0;
}
