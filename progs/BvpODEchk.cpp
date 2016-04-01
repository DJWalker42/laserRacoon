#include <cmath>
#include <string>

#include <BvpODE.h>

double model_prob_1_rhs(double x) {return 1;}
double model_prob_2_rhs(double x) {return 34. * sin(x);}

int main()
{
	std::string path = "./bin/data/";

	/* kUxx, kUx, kU, differential function, x_min, x_max */
	phys::SecondOrderODE ode_mp1(-1., 0., 0., model_prob_1_rhs, 0., 1.);

	/* lhs-type, lhs-value, rhs-type, rhs-value  */
	phys::BoundaryConditions bc_mp1(phys::DIRICHLET, 0., phys::DIRICHLET, 0.);

	/* SecondOrderODE, BoundaryConditions, number of nodes to use*/
	phys::BvpODE bvpODE_mp1(&ode_mp1, &bc_mp1, 101); 
	bvpODE_mp1.solve();
	bvpODE_mp1.plot();

	phys::SecondOrderODE ode_mp2(1., 3., -4., model_prob_2_rhs, 0., M_PI);
	phys::BoundaryConditions bc_mp2(phys::NEUMANN, -5., phys::DIRICHLET, 4.);
	phys::BvpODE bvpODE_mp2(&ode_mp2, &bc_mp2, 1001);
	bvpODE_mp2.solve();
	bvpODE_mp2.plot();

	return 0;

}
