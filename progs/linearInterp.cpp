/*
 *  Program to illustrate the limits of linear interpolation
 *  using the sinc(x) function
 */

#include "Interpolation.h"

#include <cmath> //sin()
#include <cstdio> //fprintf()

double sincSqr(double x_value);

int main(int argc, char ** argv) {

	using Linear = phys::Linear;
	using Data = phys::data;

	const int points = 100;
	const int data_points = 10;
	const double interval = 5.;

	double increment = interval / points;

	//compute the line-space x for the interval (avoiding zero)
	phys::stdVec_d x(points);
	for (int i = 0; i < points; ++i) {
		x[i] = (i + 1) * increment;
	}

	//set-up the "measured" points for the sincSqr function
	int sample {11}; //so we get 10 equidistant "measured" points
	// x indices: 0, 11, 22, 33, ...., 88, 99
	phys::stdVec_d data_x(data_points), data_y(data_points);
	for (int i = 0; i < data_points; ++i) {
		data_x[i] = x[i * sample];
		data_y[i] = sincSqr(data_x[i]);
	}

	//initialise the Linear interpolator object
	Linear linear_interp(Data(data_x, data_y));

	//interpolate the data for the given line-space x
	//How should we rewrite the either the Linear derived class or the Interpolator
	//base class so we don't have to use the class scope resolution operator to access
	//the member function we want?
	phys::stdVec_d interp_y = linear_interp.Interpolator::interpolate(x);

	//by printing to stdout we can redirect the output to a file of our choosing
	// e.g.: ./linearInterp > interp_data.log
	fprintf(stdout, "#x f(x) interp\n"); //headers for log file
	for (int i = 0; i < points; ++i)
		fprintf(stdout, "%f %f %f\n", x[i], sincSqr(x[i]), interp_y[i]);


	return 0;
}

double sincSqr(double x) {
	// limit x -> 0, sinc(x) -> 1
	return x == 0 ? 1 : sin(x ) / x ;
}
