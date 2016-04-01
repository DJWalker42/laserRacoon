#ifndef INTERPOLATION_HPP
#define INTERPOLATION_HPP

#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include <cstdio>
#include <cassert>
#include <cfloat> //DBL_EPSILON

#include "Maths.h"
#include "DynVector.h"



namespace phys{
	namespace interp{

		//Structure to coveniently handle the data to be interpolated
		struct data{
			stdVec_d x;			//!< Independent variable
			stdVec_d y;			//!< Dependent variable
			data() : x(), y(){}
			data(const stdVec_d& data_x, const stdVec_d& data_y):
				x(data_x), y(data_y){}
			void set_x_values (const stdVec_d& x_to_set){x = x_to_set;}
			void set_y_values (const stdVec_d& y_to_set){y = y_to_set;}

			void set_values(	const stdVec_d& x_to_set, 
								const stdVec_d& y_to_set	)
			{ x = x_to_set; y = y_to_set; }

			void clear(){ x.clear(); y.clear(); }
			void resize(int n){ x.resize(n); y.resize(n); }
		};

		/* Virtual base class for the different interpolation methods */
		class Interpolator{
		protected:
			Interpolator();

			Interpolator(	const stdVec_d& x_vals,
							const stdVec_d& y_vals,
							bool sort = true,
							double tolx = DBL_EPSILON );
						
			Interpolator(	const data& d,  
							bool sort = true,
							double tolx = DBL_EPSILON);
		public:
			virtual ~Interpolator(){}
		private:
			/* No Copying */
			Interpolator(const Interpolator&);
			/* No Copying */
			const Interpolator& operator=(const Interpolator&); 
		public:
			/*	Pure virtual function: interface to actually perform the interpolation
				to the point x. Derived classes implement their own version	*/
			virtual double interpolate(double x_to_interp)=0;

			/* Interpolate a set of x values; loops using the interpolate(double) interface function */
			stdVec_d interpolate(const stdVec_d& x_to_interp);

			/* Change the data you want to interpolate */
			void change_data(const data& data_to_set);

			/* Returns estimated error in the Aitken and UpDown methods only*/
			double get_error() const {return m_errorEstimate;}

		protected:
			/*	Finds the nearest data point to the x wanted.
				@returns the index, i, where x is just larger than your_data.x[i].
				If x lies outside your_data values it will return an index of zero or size()-1.
			*/
			size_t find_base_idx(double x) const;

			//checks your data is in ascending x order, if not will sort so that it is.
			void chk_sorting();		
			void sort_ascending();
		protected:
			data m_data;
			double m_errorEstimate;
			bool m_sortSwitch;
			double m_tolX;		 
		};

		/** Linear interpolation class */
		class Linear:public Interpolator{
		public:
			Linear() : Interpolator(){}
			Linear(data d, bool sort = false) : Interpolator(d, sort){}
			double interpolate(double x);
		};

		/** Lagrange interpolation class */
		class Lagrange:public Interpolator{
		public:
			Lagrange() : Interpolator(){}
			Lagrange(data d, bool sort = true) : Interpolator(d,sort){}
			double interpolate(double x);
		};

		/** Aitken method of Lagrange interpolation*/
		class Aitken:public Interpolator{
		public:
			Aitken() : Interpolator(){}
			Aitken(data d, bool sort = true) : Interpolator(d, sort){}
			double interpolate(double x);
		};

		/** UpDown method of Lagrange interpolation */
		class UpDown:public Interpolator{
		public:
			UpDown() : Interpolator(){}
			UpDown(data d, bool sort = true) : Interpolator(d, sort){}
			double interpolate(double x);
		};

		// Cubic spline interpolation
		/*	The source code for the spline interpolation was read and modified 
			from the following web address June 2014:			
				http://kluge.in-chemnitz.de/opensource/spline/ 
			Author: Tino Kluge.
		*/
		class Cubic:public Interpolator {
		public:
			Cubic():Interpolator(), m_x(), m_y(), m_a(), m_b(), m_c(), m_d(){}
			Cubic(	const data& d	):
					Interpolator(),
					m_x(), m_y(), m_a(), m_b(), m_c(), m_d()
			{
				set_data(d); // Cubic::set_data(d) -- see implementation in Spline.cpp			
			}
			void set_data(const data& d	);
			double interpolate(double x);
		private:
			stdVec_d m_x, m_y;
			stdVec_d m_a, m_b, m_c, m_d;
		};


	}
}
#endif
