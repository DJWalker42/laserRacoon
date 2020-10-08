#include <RootSearch.h>

/*
	***SHOW CASE OF ROOT SEARCHING ALGORITHMS***

	Here we use the Legendre polynomial as our search function. The derivative is required for the
	Newton-Raphson searches.

*/

/*
	This power function works for n strictly positive but what if it is called
	with n <= 0? How can we handle this? 
*/
double power(double x, int n)
{
	double retval = x;
	for(int i = 1; i<n; i++)
	{
		x *= retval;
	}
	retval = x;
	return retval;
}


double function(double x)
{
	double x2, x4, x6, x8;
	double c0 = 35, c2 = -1260, c4 = 6930, 
		c6 = -12012, c8 = 6435, d = 128;

	x2 = power(x,2);
	x4 = power(x,4);
	x6 = power(x,6);
	x8 = power(x,8);

	return ( (c8*x8 + c6*x6 + c4*x4 + c2*x2 + c0)/d );
}

double derivative(double x)
{
	double x3, x5, x7;
	double c2 = -1260, c4 = 6930, 
		c6 = -12012, c8 = 6435, d = 128;

	x3 = power(x,3);
	x5 = power(x,5);
	x7 = power(x,7);

	return ( (8*c8*x7 + 6*c6*x5 + 4*c4*x3 + 2*c2*x)/d );
}


int main ()
{
	double left = 0.3;
	double right = 0.7;

	/*
	 *  We don't have to use pointers to base class here but if you wanted to put these
	 *  objects into a single container object then this is (one way) of doing that.
	 *  How would you change this code so that it takes advantage of more modern C++
	 *  constructs? (Hint: smart pointers)
	 */

	phys::RootSearch* bisect = new phys::Bisection(&function, left, right);
	phys::RootSearch* secant = new phys::Secant(&function, left, right);
	phys::RootSearch* newton = new phys::Newton_Raphson(&function, &derivative, left);
	phys::RootSearch* hybri1 = new phys::Hybrid_B_N_R(&function, &derivative, left, right);
	phys::RootSearch* hybri2 = new phys::Hybrid_B_S(&function, left, right);
	bisect->find_root(1); //the argument requests the function to print its output
	secant->find_root(1);
	newton->find_root(1);
	hybri1->find_root(1);
	hybri2->find_root(1);
	delete bisect;
	delete secant;
	delete newton;
	delete hybri1;
	delete hybri2;

	return 0;
}
