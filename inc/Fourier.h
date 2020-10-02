#ifndef LASER_RACOON_FOURIER_HPP
#define LASER_RACOON_FOURIER_HPP

#include "DynVector.h" //contains definition of complex numbers

namespace phys{

/*	Performs a discrete Fourier transform of the data
	Use for instruction purposes only
*/
const stdVec_c DFT( const stdVec_c& );

/*	Performs a version of the fast Fourier transform of the data
	Pass as non-const reference as data has to be resized
	(padded with zeros) to be a power of two.
*/
const stdVec_c FFT( stdVec_c& );

} //namespace

#endif //header guard
