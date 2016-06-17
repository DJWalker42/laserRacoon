#include <cmath>
#include <iomanip>
#include <algorithm>

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

		if (abs_n <= zeroThresh)
		{
			ss << 0 << ".";
			for (int i = 1; i < sf; ++i) //start at 1 as "0." is already 1 s.f.
				ss << 0;
		}
		else if ((abs_n < minSci || abs_n > maxSci)) //anything less than minSci or greater than maxSci is displayed in scientific format
		{
			int d = int(ceil(log10(abs_n)));
			double order = pow(10, sf - d);
			ss << std::scientific << std::setprecision(sf - 1) << std::round(order * n) / order;
		}
		else
		{
			int d = int(ceil(log10(abs_n)));
			double order = pow(10., sf - d);

			double rn = std::round(n * order) / order;

			if (d == int(log10(rn < 0 ? -rn : rn))) ++d; //account for the rounded value being an exact power of ten

			ss << std::fixed << std::setprecision(std::max(sf - d, 0)) << rn;
		}

		return ss.str();
	}
}
