#include "Fourier.h"
#include "Maths.h" //next_power_of_2
//uncomment the following lines to use
//static const size_t log2_lut[] = { 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072};
//static const size_t log2_lut_size = sizeof(log2_lut)/sizeof(log2_lut[0]);

/*
 *  The transforms as produced by these Fourier transform functions has the following structure:
 *  for a 'data' vector of length N (== 2**k)
 *  > the "zero" frequency is located at index zero;
 *  > positive frequencies correspond to indexes 1 -> N/2 - 1;
 *  > negative frequencies correspond to indexes N/2 + 1 -> N - 1 (most negative to least negative);
 *  > and the index N/2 gives the Nyquist critical frequency (either positive or negative)
 *
 *  To test the correctness of these functions one should perform the transform on a particular
 *  'data' vector of real numbers, take the complex conjugate of the resulting transform, then
 *  perform another transform (the inverse) on that modified data. Normalising the result of the
 *  inverse transform should recover the original 'data' vector (accounting for unit round-off
 *  precision).
 *
 *  Test data to try (complex numbers i.e. (re,im)):
 *
 *  d = [ (1,0), (1,0), (1,0), (1,0), (0,0), (0,0), (0,0), (0,0)]
 *  F{d} = [(4,0), (1,-2.41421), (0,0), (1,-0.414214), (0,0), (1,0.414214), (0,0), (1,2.41421)]
 *
 *
 *  Notes: These are non-optimised algorithms, for example
 */

namespace phys {

const stdVec_c DFT(const stdVec_c &data) {
	double q;
	size_t N = data.size();

	const double pi = 4.0 * atan(1.0);
	double x = 2. * pi / double(N);

	double result_real;
	double result_imag;

	stdVec_c result;

	for (size_t i = 0; i < N; ++i) {
		result_real = 0.0;
		result_imag = 0.0;
		//complex result_tmp;
		for (size_t j = 0; j < N; ++j) {
			q = x * j * i;

			/*
			 *  this is complex number multiplication i.e.
			 *  w = complex(cos(q), sin(q))
			 *  result_tmp += data[j] * w ;
			 */

			result_real += data[j].real() * cos(q) - data[j].imag() * sin(q);
			result_imag += data[j].imag() * cos(q) + data[j].real() * sin(q);
		}
		result.push_back(complex(result_real, result_imag));
		//result.push_back( result_tmp );
	}
	return result;
}

const stdVec_c FFT(stdVec_c &data) {
	size_t N = data.size();
	size_t tN = phys::maths::next_power_of_2(N);
	N = tN;
	size_t M = 0;
	while (tN >>= 1)
		++M;

	//then pad data with zeros if necessary
	data.resize(N, complex(0., 0.));

	double pi = 4.0 * atan(1.0);

	size_t N2 = N / 2, j;

	//intial value for p has to be one.
	size_t p = 1;

	//take a copy of the data.
	stdVec_c result(data);

	//Rearrange data into bit reversed order
	for (size_t k = 1; k < N; ++k) {
		if (k < p) {
			result[p - 1] = data[k - 1];
			result[k - 1] = data[p - 1];
		}
		j = N2;
		while (j < p) {
			p -= j;
			j /= 2;
		}
		p += j;
	}

	size_t p2 = 1, p1;
	double q, u, v;
	double rtmp, itmp;

	//Apply the required additions for the transform
	for (size_t p = 0; p < M; ++p) {
		q = 0.0;
		p1 = p2;
		p2 = 2 * p1;
		for (size_t k = 0; k < p1; ++k) {
			u = cos(q);
			v = -sin(q);
			//complex w(cos(q),sin(q))
			q += pi / p1;
			for (size_t j = k; j < N; j += p2) {

				/*
				 *  complex temp = result[j + p1] * w;
				 *  result[j + p1] = result[j] - temp;
				 *  result[j] = result[j] + temp;
				 */

				rtmp = result[j + p1].real() * u - result[j + p1].imag() * v;
				itmp = result[j + p1].real() * v + result[j + p1].imag() * u;

				result[j + p1] = complex(result[j].real() - rtmp, result[j].imag() - itmp);
				result[j] = complex(result[j].real() + rtmp, result[j].imag() + itmp);
			}
		}
	}
	return result;
}

}//namespace
