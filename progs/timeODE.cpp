#include "ODE_solvers.hpp"
#include "Storage.hpp"

#include <iostream>
#include <ctime>
#include <ratio>
#include <chrono>

int main()
{
	phys::diffs::Diff_eqn* shm_eqn = new phys::diffs::SHM_eqn(9.0, 2.0);

	double initial_time = 0.0;
	double initial_pos  = 1.0;
	double initial_vel  = 0.0; 
	double step = 0.001;
	double end = 10.0;

	phys::ode::state initial_system(initial_time, initial_pos, initial_vel);  

	phys::storage::Storage container(initial_system);

	phys::ode::ODE_solver* rk4 = new phys::ode::RK4(shm_eqn, initial_system, step);

	{
		using namespace std::chrono;
		steady_clock::time_point t1 = steady_clock::now();
		while(rk4->current.x < end)
		{
			container.store(rk4->current);
			rk4->current = rk4->solve();
		}
		steady_clock::time_point t2 = steady_clock::now();
		duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
		std::cout << "Execution time: " << time_span.count() << " seconds.\n";
	}

	container.write("../Data/dummy1.txt");
	//reset system to initial values.
	rk4->set_state(initial_system);
	{
		using namespace std::chrono;
		steady_clock::time_point t1 = steady_clock::now();
		phys::storage::Storage full_result = rk4->full_solve(end);
		steady_clock::time_point t2 = steady_clock::now();
		duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
		std::cout << "Execution time: " << time_span.count() << " seconds" << std::endl;
		full_result.write("../Data/dummy2.txt");
	}
	delete rk4;
	delete shm_eqn;
	getchar();
	return 0;
}