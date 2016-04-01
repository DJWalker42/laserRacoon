#include "Maths.h"
#include <cfloat>

namespace phys{
	namespace maths{

		/** Returns a vector containing fabs(v[i]) */
		stdVec_d fabs_vector ( const stdVec_d& v )
		{ 
			stdVec_d retval(v);
			//stdVec_d::iterator it = retval.begin();
			for(auto it = retval.begin(); it != retval.end(); ++it)
				*it = fabs(*it);
			return retval;
		}

		/**	magnitude (2-norm) of two orthogonal values - avoids overflow 
			by using the maximum absolute value as a divisor.
		*/
		double mag2D(const double x, const double y)
		{
			double max, min;
			double xabs = fabs(x);
			double yabs = fabs(y);
			if( xabs > yabs )
			{
				max = xabs; 
				min = yabs;
			}
			else
			{
				max = yabs;
				min = xabs;
			}
			if(min == 0.0) return max;
			double u = min/max;
			return max * sqrt(1 + u*u);
		}

		/**	magnitude (2-norm) of three orthogonal values - avoids overflow 
			by using the maximum absolute value as a divisor.
		*/
		double mag3D(const double x, const double y, const double z)
		{
			double xabs = fabs(x);
			double yabs = fabs(y);
			double zabs = fabs(z);
			double m = max( xabs, max(yabs, zabs) );
			if (m == 0.0) return m;
			double r = m * sqrt((xabs / m) * (xabs / m) +
                    (yabs / m) * (yabs / m) +
                    (zabs / m) * (zabs / m));
			return r;
		}

		/** 2-norm of a vector of any size - avoids overflow 
			by using the maximum absolute value as a divisor.
		*/
		double magVec ( const stdVec_d& v )
		{
			if(v.empty())		return 0.0;
			if(v.size() == 1)	return fabs(v[0]);
			if(v.size() == 2)	return mag2D(v[0],v[1]);
			if(v.size() == 3)	return mag3D(v[0],v[1],v[2]);
			stdVec_d vabs = fabs_vector(v);
			double m = *std::max_element(vabs.begin(), vabs.end());
			if(m == 0.0) return m;
			double r = 0.0;
			for(unsigned i = 0; i < vabs.size(); ++i)
				r += ( (vabs[i]/m) * (vabs[i]/m) );
			return (m * sqrt(r));
		}

		/** The mean average of the values contained in the vector<double> v */
		double meanVec( const stdVec_d& v )
		{
			size_t n = v.size();
			if(n == 0) return 0.0;			
			if(n == 1) return v[0]; 
			double sum = 0.0;
			for(size_t i = 0; i < n; ++i)
				sum += v[i];
			return sum/double(n);
		}

		/* The standard deviation of values contained in a vector */
		double stdDevVec( const stdVec_d& v , double mu ) 
		{
			size_t n = v.size();
			if(n < 2) return DBL_MAX;
			double sum = 0.;
			for(size_t i = 0; i < n; ++i)
				sum += (v[i] - mu)*(v[i] - mu);
			return sqrt(sum/n);
		}

		/* Returns the mean and standard deviation of the values contained in the vector */
		std::pair<double,double> mean_stdDev( const stdVec_d& v )
		{
			double mean = meanVec(v);
			double sigma = stdDevVec(v, mean);
			return std::pair<double, double> (mean, sigma);
		}

		///** Finds the minimum and maximum values in a vector<double> v 
		//	Returns a pair of indices; 
		//	first: location of minVal in v; second: loaction of maxVal in v.
		//*/
		//std::pair<size_t,size_t> minMax(	const stdVec_d& v, 
		//									double &minVal,		
		//									double &maxVal	)
		//{		
		//	stdVec_d::const_iterator it_min = std::min_element(v.begin(), v.end());
		//	stdVec_d::const_iterator it_max = std::max_element(v.begin(), v.end());
		//	minVal = *it_min;
		//	maxVal = *it_max;
		//	std::pair<size_t,size_t> retval;
		//	retval.first = it_min - v.begin();
		//	retval.second = it_max - v.begin();
		//	return retval;
		//}

		///** Finds the minimum and maximum values in a vector<int> v 
		//	Returns a pair of indices; 
		//	first: location of minVal in v; second: location of maxVal in v.
		//*/
		//std::pair<size_t,size_t> minMax(	const stdVec_i& v, 
		//									int &minVal, 
		//									int &maxVal	)
		//{		
		//	stdVec_i::const_iterator it_min = std::min_element(v.begin(), v.end());
		//	stdVec_i::const_iterator it_max = std::max_element(v.begin(), v.end());
		//	minVal = *it_min;
		//	maxVal = *it_max;
		//	std::pair<size_t,size_t> retval;
		//	retval.first = it_min - v.begin();
		//	retval.second = it_max - v.begin();
		//	return retval;
		//}

		/**	Rounds the value to the nearest whole number: x.5 rounds to x+1. 
		*/
		double round( double val )
		{
			double flr = floor(val);
			double cil = ceil(val);
			
			double diff_flr = fabs(val - flr);
			double diff_cil = fabs(val - cil);

			if(diff_flr < diff_cil)
				return flr;
			else
				return cil;
		}

		int factorial(int n)
		{
			int result = 1;
			for(int i=1; i<=n; i++)
				result *= i;
			return result;
		} /* compute n! using a loop */

		size_t next_power_of_2(size_t n)
		{
			//find next largest power of two
			//FIX ME: what happens if n = zero is passed in? How about n = one? 
			n--;
			n |= n >> 1;
			n |= n >> 2;
			n |= n >> 4;
			n |= n >> 8;
			n |= n >> 16;
			n++;

			return n;
		}

	}
}
