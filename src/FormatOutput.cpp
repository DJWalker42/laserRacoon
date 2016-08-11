#include <sstream>
#include <cmath> //for pow, log10, and round functions
#include <iomanip> //for i/o manipulators e.g fixed
#include <algorithm> //for max function

#include "FormatOutput.h"

namespace phys{

	std::string format(	double n, 
						int sf, 
						double minSci,
						double maxSci,
						double zeroThresh)
	{
		std::stringstream ss;

		double abs_n = fabs(n);

		if (abs_n <= zeroThresh){
			ss << 0 << ".";
			for (int i = 1; i < sf; ++i) //start at 1 as "0." is already 1 s.f.
				ss << 0;

			return ss.str();
		}
		
		int d = int(ceil(log10(abs_n)));
		double order = pow(10., sf - d);
		double rn = std::round(n * order) / order;
		
		//anything less than minSci or greater than maxSci is displayed in scientific format
		if ((abs_n < minSci || abs_n > maxSci)) { 
			ss << std::scientific << std::setprecision(sf - 1) << rn;
		}
		else {
			ss << std::fixed << std::setprecision(std::max(sf - d, 0)) << rn;
		}

		return ss.str();
	}
}
