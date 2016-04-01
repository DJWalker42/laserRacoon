#include <cassert>
#include <cfloat>

#include "Interpolation.h"
#include "LinearSolvers.h"

/*	The source code for the spline interpolation, including the Band Matrix definition,
	was read and modified from the following web address June 2014:			
		http://kluge.in-chemnitz.de/opensource/spline/ 
	Author: Tino Kluge.
*/

namespace phys{
	namespace interp{

		// spline implementation
		// -----------------------
		void Cubic::set_data(	const data& d	) {
			
			/*	Override sort switch to ensure data is sort ascended */
			m_sortSwitch = true;
			Interpolator::change_data(d);

			m_x = m_data.x;
			m_y = m_data.y;
			size_t n = m_x.size();

			phys::Band_matrix A(n, 1, 1);
			std::vector<double>  rhs(n);
			for (size_t i = 1; i < n - 1; i++) {
				A(i, i - 1) = 1.0 / 3.0*(m_x[i] - m_x[i - 1]);
				A(i, i) = 2.0 / 3.0*(m_x[i + 1] - m_x[i - 1]);
				A(i, i + 1) = 1.0 / 3.0*(m_x[i + 1] - m_x[i]);
				rhs[i] = (m_y[i + 1] - m_y[i]) / (m_x[i + 1] - m_x[i]) - (m_y[i] - m_y[i - 1]) / (m_x[i] - m_x[i - 1]);
			}
			// boundary conditions, zero curvature b[0]=b[n-1]=0
			A(0, 0) = 2.0;
			A(0, 1) = 0.0;
			rhs[0] = 0.0;
			A(n - 1, n - 1) = 2.0;
			A(n - 1, n - 2) = 0.0;
			rhs[n - 1] = 0.0;

			// solve the equation system to obtain the parameters b[]
			m_b = A.lu_solve(rhs);

			// calculate parameters a[] and c[] based on b[]
			m_a.resize(n);
			m_c.resize(n);
			for (size_t i = 0; i < n - 1; i++) {
				m_a[i] = 1.0 / 3.0*(m_b[i + 1] - m_b[i]) / (m_x[i + 1] - m_x[i]);
				m_c[i] = (m_y[i + 1] - m_y[i]) / (m_x[i + 1] - m_x[i])
					- 1.0 / 3.0*(2.0*m_b[i] + m_b[i + 1])*(m_x[i + 1] - m_x[i]);
			}

			// for the right boundary we define
			// f_{n-1}(x) = b*(x-x_{n-1})^2 + c*(x-x_{n-1}) + y_{n-1}
			double h = m_x[n - 1] - m_x[n - 2];
			// m_b[n-1] is determined by the boundary condition
			m_a[n - 1] = 0.0;
			m_c[n - 1] = 3.0*m_a[n - 2] * h*h + 2.0*m_b[n - 2] * h + m_c[n - 2];   // = f'_{n-2}(x_{n-1})
		}

		double Cubic::interpolate (double x){			
			// find the closest point m_x[idx] < x, idx = 0 even if x < m_x[0]
			size_t i0 = find_base_idx(x);

			/*	quick return if x within tolerance of your data x value */
			if (fabs(x - m_data.x[i0]) < m_tolX) return m_data.y[i0];

			size_t n = m_x.size();
			double h = x - m_x[i0];
			double interpol;
			if (x < m_x[0]) {
				// extrapolation to the left
				interpol = ((m_b[0])*h + m_c[0])*h + m_y[0];
			}
			else if (x > m_x[n - 1]) {
				// extrapolation to the right
				interpol = ((m_b[n - 1])*h + m_c[n - 1])*h + m_y[n - 1];
			}
			else {
				// interpolation
				interpol = ((m_a[i0] * h + m_b[i0])*h + m_c[i0])*h + m_y[i0];
			}
			return interpol;
		}

	}
}
