#include <ODESolvers.h>
#include <Visualise.h>

/*
	Testing some different ODE solvers on a damped spring. 

	Take note that the Stormer-Verlet method is designed to work with 
	differential equations containing only second derivatives of the 
	dependent variable on the right-hand-side, i.e. the damping term 
	is irrelevant. It will compute a solution with no damping i.e.
	SHM for the given spring constant.

	Note also that this means the first derivative storage vector will
	contain only zeros for the Stormer-Verlet solution rendering a
	phase-space plot meaningless. A message is printed to the console
	to that affect. 

*/

int main()
{
	const int num_of_solvers = 7;
	std::string root = "./data/";

	phys::Diff_eqn* shm_func = new phys::SHM(9.0, 1.0);//spring const., drag coeff.
	
	double step = 0.01;
	double initial_time = 0.0;
	double initial_position = 1.0;
	double initial_velocity = 0.0; 

	phys::state initial_system(initial_time, initial_position, initial_velocity);

	std::vector<phys::ODEStorage> container(num_of_solvers);
	std::vector<phys::ODESolver*> solvers(num_of_solvers);
	std::vector<std::string> filenames(num_of_solvers);
	std::vector<std::string> titles(num_of_solvers);

	solvers[0] = new phys::Euler(shm_func, initial_system, step);
	solvers[1] = new phys::Imp_Euler(shm_func, initial_system, step);
	solvers[2] = new phys::Mod_Euler(shm_func, initial_system, step);
	solvers[3] = new phys::Velocity_Verlet(shm_func, initial_system, step);
	solvers[4] = new phys::Stormer_Verlet(shm_func, initial_system, step);
	solvers[5] = new phys::RK4(shm_func, initial_system, step);
	solvers[6] = new phys::Leapfrog(shm_func, initial_system, step);

	filenames[0] = root + "euler_shm.txt";
	filenames[1] = root + "imp_euler_shm.txt";
	filenames[2] = root + "mod_euler_shm.txt";
	filenames[3] = root + "verlet_shm.txt";
	filenames[4] = root + "sto_verlet_shm.txt";
	filenames[5] = root + "rk4_shm.txt";
	filenames[6] = root + "leapfrog_shm.txt";

	titles[0] = "Euler";
	titles[1] = "Improved Euler";
	titles[2] = "Modified Euler";
	titles[3] = "Velocity Verlet";
	titles[4] = "Stormer Verlet";
	titles[5] = "Runge Kutta";
	titles[6] = "Leapfrog";

	phys::Viewer viewer;
	for(unsigned i = 0; i < container.size(); i++)
	{
		container[i] = solvers[i]->fullSolve(10.);
		//container[i].write(filenames[i]);
		viewer.set_plot_name(titles[i]);
		try{
			viewer.plot(container[i], phys::Viewer::PHASE); //use Viewer::DEPEN to see y(t) plots instead
		} catch(const std::exception& e) {
			if (i == 4){
				std::cout << "Cannot plot phase-space for the Stormer-Verlet solution as it contains no velocity data\n";
			} else {
				std::cout << e.what() << std::endl;
			}
		}
		viewer.clear();
	}

	//to test ODEStorage::read member function elsewhere
	container[0].write("./euler_shm.log", true);


	return 0;
}
