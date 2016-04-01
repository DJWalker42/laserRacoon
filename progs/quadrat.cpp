#include <Quadrature.h>
#include <iomanip>

double f (double x)
{
	return 1./(x*x + 1.);
}


int main()
{
	double a = 0., b = 6.; 

	phys::quad::Quadrature* p_mido = new phys::quad::Mid_ordinate;
	phys::quad::Quadrature* p_trap = new phys::quad::Trapeziod;
	phys::quad::Quadrature* p_simp = new phys::quad::Simpson;
	phys::quad::Quadrature* p_bool = new phys::quad::Boole(true);

	std::cout << "Gauss-Legendre\n";
	double sol = 0.;
	phys::quad::Legendre legen(2);
	for (int i = 3; i < 33; ++i)
	{
		sol = legen.integrate(f, a, b);
		legen.set_points(i);
		std::cout << i - 1 << " : " << std::setprecision(16) << sol << "\n";
	}
	std::cout << " --\n";	

	std::cout << "Strips\t" << "Sol\t" << "Error\t" << "est. Error\n";
	std::cout << "Mid-Ordinate\n";
	for(int i = 1; i<11; ++i)
	{
		double r = p_mido->integrate(f, a, b, 2 * i);
		std::cout << 2*i << "\t" << r <<"\t" 
			<< sol - r << "\t" << p_mido->get_error() << "\n";  
	}
	std::cout << " --\n";

	std::cout << "Trapezoid\n";
	for(int i = 1; i<11; ++i)
	{
		double r = p_trap->integrate(f, a, b, 2 * i);
		std::cout << 2 * i << "\t" << r << "\t"
			<< sol - r << "\t" << p_trap->get_error() << "\n";
	}
	std::cout << " --\n";

	std::cout << "Simpson\n";
	for(int i = 1; i<11; ++i)
	{
		double r = p_simp->integrate(f, a, b, 2 * i);
		std::cout << 2 * i << "\t" << r << "\t"
			<< sol - r << "\t" << p_simp->get_error() << "\n";
	}
	std::cout << " --\n";

	std::cout << "Boole\n";
	for(int i = 1; i<11; ++i)
	{
		double r = p_bool->integrate(f, a, b, 2 * i);
		std::cout << 4 * i << "\t" << r << "\t"
			<< sol - r << "\t" << p_bool->get_error() << "\n";
	}
	std::cout << " --\n";

	p_simp->reset_func_calls(); 
	phys::quad::Adaptive adapt(p_simp, 1.e-6);
	double r = adapt.integrate(f, a, b );
	std::cout << "Adaptive: " << r << " Error = " << sol - r << std::endl;
	std::cout << adapt.get_f_calls() << "\n"; 

	phys::quad::Romberg romb(1.e-6);
	double r1 = romb.integrate(f, a, b);

	std::cout << "Romberg: " << r1 << " Error = " << sol - r1 << " Err = " << romb.get_error() << std::endl;

	std::cout << romb.get_func_calls() << std::endl;

	return 0;
}
