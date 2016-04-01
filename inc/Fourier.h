#ifndef FOURIER_HPP
#define FOURIER_HPP

#include "Complex.h"
#include "DynVector.h"

namespace phys{
	namespace fourier{

		/*	Performs a discrete fourier transform of the data 
			Use for instruction purposes only
		*/
		const stdVec_c DFT( const stdVec_c& ); 

		/*	Performs a version of the fast fourier transform of the data
			Pass as non-const reference as data has to be resized 
			(padded with zeros) to be a power of two.
		*/
		const stdVec_c FFT( stdVec_c& );

	}

}


#endif
