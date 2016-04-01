#include <Quadrature.h>
#include <Visualise.h>

double pi = 4.0*atan(1.0);

double f1(double x){
	return cos(pi*x*x/2);
}

double f2(double x){
	return sin(pi*x*x/2);
}


int fresnel_guass()
{
	phys::quad::Legendre lege(20);

	//set lower limit to zero, increment upper limit to a max of 5.0
	double b = 0.0;
	double T1, T2;
	double delta = 0.05;
	phys::stdVec_d v, I;

	while (b < 5.0)
	{
		b += delta;
		T1 = lege.integrate(f1, 0.0, b) + 0.5;
		T2 = lege.integrate(f2, 0.0, b) + 0.5;
		I.push_back((T1*T1 + T2*T2) / 2);
		v.push_back(b);
	}

	

	phys::visual::Viewer viewer;

	viewer.withLines();

	viewer.set_plot_name("Guass");

	viewer.plot(v, I);

	return lege.get_func_calls(); 

}

int fresnel_romberg()
{

	phys::quad::Romberg romb;

	//set lower limit to zero, adjust upper limit to a max of 5.0
	double b = 0.0;
	double T1, T2;
	double delta = 0.05;
	std::vector<double> v, I;

	while (b <= 5.0)
	{
		b += delta;
		T1 = romb.integrate(f1, 0.0, b) + 0.5;
		T2 = romb.integrate(f2, 0.0, b) + 0.5;
		I.push_back((T1*T1 + T2*T2) / 2);
		v.push_back(b);
	}

	phys::visual::Viewer viewer;

	viewer.set_plot_name("Romberg");
	viewer.withLines();

	viewer.plot(v, I);

	return romb.get_func_calls();
}



int main ()
{
	int rombe_num = fresnel_romberg();
	int guass_num = fresnel_guass();

	std::cout << "Romberg function calls: " << rombe_num << "\n"; 
	std::cout << "Legendre function calls: " << guass_num << "\n";
	
	return 0;
}
