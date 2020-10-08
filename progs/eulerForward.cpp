#include "ODESolvers.h"

double phys::User_eqn::differential_function(double x, const stdVec_d& y, int N, int c ) {
	return -x * y[0];
}

double solution(double x) {
	return exp(-x * x / 2.);
}

int main (int argc, char ** argv) {

	using Euler = phys::Euler;
	using Differential = phys::User_eqn;
	using ODEStorage = phys::ODEStorage;


	Differential df(1); //differential is of order 1
	double x_start = 0.0;
	double x_end = 2.0;
	double y0 = 1.0;
	phys::state initial_state(x_start, y0);

	double soln = solution(x_end);
	std::cout << "Soln: " << soln << "\n\n";

	for (int n = 2; n < 8; n++) {

		//half the step size from previous iteration
		double step = (x_end - x_start) / pow(2.,n);

		Euler euler_solver(&df, initial_state, step);

		ODEStorage storage = euler_solver.fullSolve(x_end);

		//get the computed value of y at x_end
		double euler_soln = storage.get_dependent().back();

		std::cout << "Step: " << step;
		std::cout << "\tEuler: " << euler_soln;
		std::cout << "\terror: " << fabs(soln - euler_soln);
		std::cout << "\n";
	}

	std::cout << std::endl;

	return 0;
}
