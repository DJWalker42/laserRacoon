#include <Quadrature.h>
#include <Visualise.h>
#include <cmath>

#ifndef PI
	#define PI 3.1415926535897932384626433832795028841971693993751 //slight overkill on pi precision
#endif

/*	List of standard integrals and their solutions

	Integrand				Solution

	sqrt(x) * exp(-x);		sqrt(PI)/2.;
	exp(-a*x*x);			sqrt(PI/a)/2.;		// a > 0
	x*x*exp(-a*x*x);		sqrt(PI/a/a/a)/4.:	// a > 0
	x*x*x*exp(-a*x*x);		1/2./a/a;			// a > 0
	x/(exp(x) - 1);			PI*PI / 6.;
	x*x/(exp(x) - 1);		2.40;				//approx.
	x*x*x/(exp(x) - 1);		PI^4 / 15.;

	x ? sin(x)/x : 1;			PI/2.;
	x ? sin(x)*sin(x)/x/x : 1;	PI/2.;

*/

double I0(double x) { return sqrt(x) * exp(-x); }
double I1(double x) { return exp(-x*x); } 
double I2(double x) { return x*x*exp(-x*x); }
double I3(double x) { return x*x*x*exp(-x*x); }
double I4(double x) { return x / (exp(x) - 1); }
double I5(double x) { return x*x / (exp(x) - 1); }
double I6(double x) { return x*x*x / (exp(x) - 1); }
double I7(double x)	{ return x ? sin(x) / x : 1;  }
double I8(double x) { return x ? sin(x)*sin(x) / x / x : 1; }

double S0 = sqrt(PI) / 2.;
double S1 = sqrt(PI) / 2.;
double S2 = sqrt(PI) / 4.;
double S3 = 1 / 2.;
double S4 = PI*PI / 6.;
double S5 = 2.4;
double S6 = PI*PI*PI*PI / 15.;
double S7 = PI / 2.;
double S8 = PI / 2.;

typedef std::pair<double(*)(double), double> integral;

std::string functions[] = { "sqrt(x) * exp(-x)", 
							"exp(-x*x)", 
							"x*x*exp(-x*x)", 
							"x*x*x*exp(-x*x)", 
							"x / (exp(x) - 1)", 
							"x*x / (exp(x) - 1)", 
							"x*x*x / (exp(x) - 1)", 
							"sin(x) / x", 
							"sin(x)*sin(x) / x / x"
};


int main(){
	std::vector<integral> test(9);

	test[0].first = I0; test[0].second = S0;
	test[1].first = I1; test[1].second = S1;
	test[2].first = I2; test[2].second = S2;
	test[3].first = I3; test[3].second = S3;
	test[4].first = I4; test[4].second = S4;
	test[5].first = I5; test[5].second = S5;
	test[6].first = I6; test[6].second = S6;
	test[7].first = I7; test[7].second = S7;
	test[8].first = I8; test[8].second = S8;

	phys::quad::Gauss* p_lag = new phys::quad::Laguerre; 

	phys::stdVec_d numKnots;
	phys::stdVec_d value;

	for(size_t i = 0; i < 9; ++i){
		for(size_t j = 1; j < 101; ++j){
			p_lag->set_points(j);
			double result = p_lag->integrate(test[i].first);
			numKnots.push_back(j);
			value.push_back(result);
		}
		phys::visual::Viewer viewer;
		viewer.set_plot_name(functions[i]);	
		viewer.plot(numKnots, value);
		viewer.draw_line(0., test[i].second); 
	}

	return 0;
}
