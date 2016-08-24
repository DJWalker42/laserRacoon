#include <cfloat>

#include "Interpolation.h"

namespace phys{
	namespace interp{

		/********************************************************
		*	Constructors
		********************************************************/
		Interpolator::Interpolator() :	m_data(),
										m_errorEstimate(0.0),
										m_sortSwitch(true),
										m_tolX(DBL_EPSILON)
		{ /*no data to check*/}

		Interpolator::Interpolator(	const stdVec_d& x_vals,
									const stdVec_d& y_vals,
									bool sort,
									double tolx) :
									m_data(x_vals, y_vals),
									m_errorEstimate(0.0),
									m_sortSwitch(sort),
									m_tolX(tolx)
		{
			assert(m_data.x.size() == m_data.y.size());
			if (m_sortSwitch)
				chk_sorting();
		}

		Interpolator::Interpolator(	const data& d,
									bool sort,
									double tolx) :
									m_data(d),
									m_errorEstimate(0.0),
									m_sortSwitch(sort),
									m_tolX(tolx)
		{
			assert(m_data.x.size() == m_data.y.size());
			if (m_sortSwitch)
				chk_sorting();
		}

		/***************************************************************
		*	Interface functions
		***************************************************************/
		void Interpolator::chk_sorting()
		{
			size_t i = 0;
			size_t max = m_data.x.size();
			while (++i < max && m_data.x[i] > m_data.x[i - 1])
				; //no body
			if (i < max) sort_ascending();
		}

		void Interpolator::change_data(const data& data_to_set)
		{
			m_data = data_to_set;
			if (m_sortSwitch)
				chk_sorting();
		}

		stdVec_d Interpolator::interpolate(const stdVec_d& x)
		{
			stdVec_d retval;
			for(size_t i = 0; i < x.size(); ++i)
			{
				retval.push_back(this->interpolate(x[i]));
			}
			return retval;
		}

		void Interpolator::sort_ascending()
		{
			data temp(m_data);
			//this sorts m_data x values.
			std::sort(m_data.x.begin(), m_data.x.end());
			//need to apply to y values.
			for (size_t i = 0; i < temp.x.size(); ++i)
			{
				size_t j = 0;
				while( j < m_data.x.size() && temp.x[i] != m_data.x[j] )
					++j; //put update in the body to ensure we know which condition breaks the loop

				if( j < m_data.x.size() )
					temp.y[j] = m_data.y[i];
				//else no match - which is impossible as temp is a copy of m_data
			}
			//assign sorted data back to m_data.
			m_data.y = temp.y;
		}

		size_t Interpolator::find_base_idx(double x) const
		{
			// find the closest point m_data.x[idx] < x, idx = 0 even if x < m_x[0]
			stdVec_d::const_iterator it;
			it = std::lower_bound(m_data.x.begin(), m_data.x.end(), x);
			return std::max(static_cast<int>(it - m_data.x.begin() - 1), 0);
		}


		/********************************************************************************
		*	Derived Class implementations
		********************************************************************************/
		double Linear::interpolate(double x)
		{
			size_t i0 = find_base_idx(x);
			/*	quick return if x equal or near to your data x value */
			if (fabs(x - m_data.x[i0]) < m_tolX) return m_data.y[i0];

			/*	catch possible right hand extrapolation outcome */
			if(i0 == m_data.x.size() - 1) --i0;

			/*	for left hand extrapolation outcome i0 will be zero which is not a problem */
			double c1 = (x - m_data.x[i0+1])/(m_data.x[i0] - m_data.x[i0+1]);
			double c2 = (x - m_data.x[i0])/(m_data.x[i0+1] - m_data.x[i0]);
			return (c1*m_data.y[i0] + c2*m_data.y[i0+1]);
		}

		double Lagrange::interpolate(double x)
		{
			size_t i0 = find_base_idx(x);
			/*	quick return if x equal or near to your data x value */
			if (fabs(x - m_data.x[i0]) < m_tolX) return m_data.y[i0];
			size_t m = m_data.x.size();
			double retval = 0.0;
			stdVec_d lambda(m, 1.0);
			for(size_t k = 0; k < m; ++k)
			{
				for(size_t l = 0; l < m; ++l)
				{
					if(k != l) lambda[k] *=
						(x - m_data.x[l])/(m_data.x[k] - m_data.x[l]);
				}

				retval += lambda[k]*m_data.y[k];
			}
			return retval;
		}

		double Aitken::interpolate(double x)
		{
			size_t m = m_data.x.size();
			size_t i0 = find_base_idx(x);
			/*	quick return if x equal or near to your data x value */
			if (fabs(x - m_data.x[i0]) < m_tolX) return m_data.y[i0];

			double x1, x2, f1 = 0.0, f2 = 0.0;
			data temp(m_data);
			for(size_t i = 1; i < m; ++i)
			{
				for(size_t j = 0; j < m - i; ++j)
				{
					x1 = temp.x[j];
					x2 = temp.x[j+i];
					f1 = temp.y[j];
					f2 = temp.y[j+1];

					temp.y[j] = f2*(x - x1)/(x2 - x1) + f1*(x - x2)/(x1 - x2);
				}
			}
			double retval = temp.y[0];
			m_errorEstimate = (fabs(retval - f1)+fabs(retval - f2))/2.0;
			return retval;
		}

		double UpDown::interpolate(double x)
		{
			const size_t n_max = 21; //i.e. max polynomial order of 20
			size_t m = m_data.x.size();
			assert(m <= n_max);

			size_t i0 = find_base_idx(x);
			/*	quick return if x equal or near to your data x value */
			if (fabs(x - m_data.x[i0]) < m_tolX) return m_data.y[i0];

			double dx = fabs(x - m_data.x[i0]);
			double dxt = (i0 < m_data.x.size() - 1)? fabs(x - m_data.x[i0+1]) : dx;
			if(dxt < dx){++i0;}
			size_t j0 = i0;

			// set up correction matrices and explictly assign to zero.
			double dp[n_max][n_max];
			double dm[n_max][n_max];

			for(int i = 0; i < n_max; ++i){
				for(int j = 0; j < n_max; ++j){
					dp[i][j] = 0.0;
					dm[i][j] = 0.0;
				}
			}

			//Compute the correction matrices
			for(size_t i = 0; i < m; ++i)
			{
				dp[i][i] = m_data.y[i];
				dm[i][i] = m_data.y[i];
			}
			for(size_t i = 1; i < m; ++i)
			{
				for(size_t j = 0; j < m - i; ++j)
				{
					size_t k = i + j;
					dx = (dp[j][k-1] - dm[j+1][k])/(m_data.x[k] - m_data.x[j]);
					dp[j][k] = dx*(m_data.x[k] - x);
					dm[j][k] = dx*(m_data.x[j] - x);
				}
			}

			//update the approximation
			double retval = m_data.y[i0];
			size_t it = (x < m_data.x[i0]) ? 1 : 0;
			for(size_t i = 1; i < m; ++i)
			{
				if( it == 1 || j0 == m - 1)
				{
					i0 -= 1;
					m_errorEstimate = dp[i0][j0];
					retval += m_errorEstimate;
					it = (j0 == m - 1) ? 1 : 0;
				}else if( it == 0 || i0 == 0 )
				{
					j0 += 1;
					m_errorEstimate = dm[i0][j0];
					retval += m_errorEstimate;
					it = (i0 == 0) ? 0 : 1;
				}
			}
			m_errorEstimate = fabs(m_errorEstimate);
			return retval;
		}

		/*************************************************************************
		*	Cubic Spline Implementation in "Spline.cpp"
		*************************************************************************/


	}
}
