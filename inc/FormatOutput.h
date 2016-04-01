#ifndef FORMATOUTPUT_HPP
#define FORMATOUTPUT_HPP

#include <sstream>

#ifndef DBL_EPSILON
#define DBL_EPSILON 2.2204460492503131E-16
#endif

namespace phys{

	//return a string of the number supplied to the number of significant figures requested
	//optional - specify the minimum absolute number below which the output is in scientific notation
	//optional - specify the maximum absolute number above which the output is in scientific notation
	//Numbers less than the threshToZero argument will be returned as 0.0... up to the number of sig. figs. requested.
	std::string format(	double number, 
						int significantFigs, 
						double minNumScientific = 1.e-2,
						double maxNumScientific = 1.e2,
						double threshToZero = DBL_EPSILON);
}

#endif
