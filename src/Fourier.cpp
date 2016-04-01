#include "Fourier.h"
#include "Maths.h"

static const size_t log2_lut[] = { 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072}; 
static const size_t log2_lut_size = sizeof(log2_lut)/sizeof(log2_lut[0]);


namespace phys{
	namespace fourier{

		const stdVec_c DFT( const stdVec_c& data )
		{
			double q;
			size_t N = data.size();

			const double pi = 4.0*atan(1.0);
			double x = 2*pi/double(N);

			double result_real;
			double result_imag;

			stdVec_c result;

			for(size_t i = 0; i < N; ++i)
			{
				result_real = 0.0;
				result_imag = 0.0;
				for(size_t j = 0; j < N; ++j)
				{
					q = x*j*i;
					// complex is a point2 type structured as ( real, imag ).
					result_real += data[j][0]*cos(q) + data[j][1]*sin(q);
					result_imag += data[j][1]*cos(q) - data[j][0]*sin(q);
				}
				result.push_back( complex(result_real, result_imag) );
			}
			return result;
		}

		const stdVec_c FFT( stdVec_c& data )
		{
			size_t N = data.size();
			size_t tN = phys::maths::next_power_of_2(N);
			N = tN;
			size_t M = 0;
			while(tN >>= 1)
				++M;

			//then pad data with zeros if nescessary
			data.resize(N, complex(0.,0.));

			double pi = 4.0*atan(1.0);
	
			size_t N2 = N/2 , j;

			//intial value for p has to be one.
			size_t p = 1;

			//take a copy of the data.
			stdVec_c result(data);

			//Rearrange data into bit reversed order
			for(size_t k = 1; k < N; ++k)
			{
				if(k < p )
				{
					result[p-1] = data[k-1];
					result[k-1] = data[p-1];
				}
				j = N2;
				while(j < p)
				{
					p -= j;
					j /= 2;
				}
				p += j;
			}

			size_t p2 = 1, p1;
			double q, u, v;
			double rtmp, itmp;

			//Apply the required addtions for the transform
			for(size_t p = 0; p < M; ++p)
			{
				q = 0.0;
				p1 = p2;
				p2 = 2*p1;
				for(size_t k = 0; k < p1; ++k)
				{
					u = cos(q);
					v = -sin(q);
					q += pi/p1;
					for(size_t j = k; j < N; j += p2)
					{
						rtmp = result[j + p1][0]*u - result[j + p1][1]*v;
						itmp = result[j + p1][0]*v + result[j + p1][1]*u;

						result[j + p1] = phys::complex( result[j][0] - rtmp, result[j][1] - itmp );
						result[j] = phys::complex( result[j][0] + rtmp, result[j][1] + itmp ); 
					}
				}
			}
			return result;
		}

	}

}
