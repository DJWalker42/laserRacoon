#include <ODESolvers.h>
#include <Visualise.h>

#include <iomanip>

/*
	***SIMULATION OF A DAMPED DRIVEN OSCILLATOR***
	
	This program shows how the response of a driven oscillator is affected by damping.

	How might you obtain the Q-factor from the data computed?

*/

//set global driving frequency to use in the the drive function
//and which can be modified in the main function
const double omega_start = 1.;
double omega = omega_start;

double drive_func(double t)
{
	return cos(omega*t);
}

int main()
{
	double k = 9.; //natural frequency w0 = 3 Hz; mass == 1.
	double D = .5;

	phys::diffs::SHM shm_eq(k, D, drive_func);
	phys::ode::state initial(0., 1., 0.);
	double step = .01;
	phys::ode::RK4 rk4_solver(&shm_eq, initial, step);

	double omega_step = .05, omega_max = 5.;
	double drag_step = .5, drag_max = 3.;

	size_t omega_N = size_t( (omega_max - omega_start)/omega_step );
	size_t drag_N = size_t( (drag_max - D)/drag_step );

	phys::stdVec_d omega_vals(omega_N);
	std::vector<phys::stdVec_d> amp_vals(drag_N);

	for(size_t j = 0; j < omega_N; ++j)
	{
		omega_vals[j] = omega;
		omega += omega_step;
	}

	std::vector<std::string> key_names(drag_N);

	phys::visual::Viewer viewer;
	viewer.set_pause(30);
	viewer.set_y_range(-1., 1.);

	for(size_t i = 0; i < drag_N; ++i) 
	{
		omega = omega_start;
		for(size_t j = 0; j < omega_N; ++j)
		{
			//solve up to 50 seconds to allow solution to stabilise
			phys::storage::ODEStorage data = rk4_solver.fullSolve(50.); 
			viewer.plot(data);
			viewer.clear();
			//find the max amplitude in the stable region - assume we reached stability by half way point.
			phys::stdVec_d amplitude = data.get_dependent(); 
			phys::stdVec_d sub_amp(amplitude.begin() + amplitude.size()/2, amplitude.end()); 
			double maxVal, minVal;
			phys::maths::minMax(sub_amp, minVal, maxVal);
			amp_vals[i].push_back(maxVal);
			//increase driving frequency		
			omega += omega_step;			 
			//reset solver to initial state
			rk4_solver.set_system(initial);			
		}
		//increase the drag
		D += drag_step;
		shm_eq.set_drag_coeff(D);
		std::ostringstream convert;
		convert << std::setprecision(2) << std::fixed << D; 
		key_names[i] = convert.str();
	}
	viewer.set_pause(0);
	viewer.set_y_range(0., .7);
	viewer.set_key_name( key_names );
	viewer.plot(omega_vals, amp_vals);
    //to save an image of the plot uncomment line below and give a directory location and filename (.png, .jpg, .bmp).
    //viewer.save(/*supply a location + filename*/);

	return 0;
}
