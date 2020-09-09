#ifndef PHYSMATHS_HPP
#define PHYSMATHS_HPP

#include <cmath>
#include <algorithm>

#include "DynVector.h"
#include "Complex.h"

#ifdef USING_PHYS_MATHS_DEF
	#define PI 3.1415926535897932384626433832795028841971693993751 
#endif

namespace phys{
	namespace maths{

		/** Returns a vector containing fabs(v[i]) */
		stdVec_d fabs_vector ( const stdVec_d& v );

		template<typename T>
		T max( const T a, const T b ){ return a > b ? a : b; }

		template<typename T>
		T min(const T a, const T b){ return a < b ? a : b; }

		/**	@brief computes magnitude of a (cartesian) 2-D vector 
			@param x double value x-component
			@param y double value y-component
			@returns magnitude double value.
		*/
		double mag2D(const double x, const double y);

		/**	@brief computes magnitude of a (cartesian) 3-D vector 
			@param x double value x-component
			@param y double value y-component
			@param z double value z-component
			@returns magnitude double value.
		*/
		double mag3D(const double x, const double y, const double z);

		/**	@brief computes magnitude of a general vector of doubles
			@param v vector of doubles
			@returns magnitude double value.
		*/
		double magVec ( const stdVec_d& v );

		/** @brief computes the arithmetic mean of the values contained in a vector of doubles
		*/
		double meanVec( const stdVec_d& v); 

		double stdDevVec( const stdVec_d& v , double mu );

		std::pair<double,double> mean_stdDev( const stdVec_d& v ); 


		template<typename T>
		std::pair<size_t, size_t> minMax(	const std::vector<T>& values,
											T& minVal,
											T& maxVal )
		{
			typename std::vector<T>::const_iterator it_min = std::min_element(values.begin(), values.end());
			typename std::vector<T>::const_iterator it_max = std::max_element(values.begin(), values.end());
			minVal = *it_min;
			maxVal = *it_max;
			std::pair<size_t, size_t> retval;
			retval.first = it_min - values.begin();
			retval.second = it_max - values.begin();
			return retval;
		}


		///**	@brief locates minimum and maximum values contained in a vector of doubles
		//	@param minVal reference to the minimum value found (double)
		//	@param maxVal reference to the maximum value found (double)
		//	@returns a pair of indices: first is location of minVal, second is location of maxVal.
		//*/
		//std::pair<size_t,size_t> minMax(	const stdVec_d& values, 
		//									double &minVal, 
		//									double &maxVal );
		///**	@brief locates minimum and maximum values contained in a vector of integers
		//	@param minVal reference to the minimum value found (int)
		//	@param maxVal reference to the maximum value found (int)
		//	@returns a pair of indices: first is location of minVal, second is location of maxVal.
		//*/
		//std::pair<size_t,size_t> minMax(	const stdVec_i& values, 
		//									int &minVal, 
		//									int &maxVal);
		/** @brief rounds the value to the nearest whole digit with the policy that X.5 rounds up to X+1
		*/
		double round ( double val );

		/** Computes the factorial (!) of n using a loop*/
		int factorial(int n);

		template<class T>
		const T pow(const T&x, int n)
		{
			T result = 1;
			T powerOfX = x;
			while(n){
				if(n % 2) result *= powerOfX;
				powerOfX *= powerOfX;
				n /= 2;
			}
			return result;
		} // compute a power using a loop

		template<class T>
		const T power(const T&x, int n)
		{
			return n ? (n%2 ? x * power(x * x,n/2) : power(x * x,n/2)) : 1;
		} // compute a power recursively

		//@param	n	the number from which we want to find the next largest power of two			
		//@return		The next largest power of two
		size_t next_power_of_2(size_t n);

	}
}



#endif
