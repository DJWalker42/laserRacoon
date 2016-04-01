#include <physMaths.hpp>
#include <iostream>


namespace phys{

	class interval{
	public:
		interval();
		interval(double s, double e);

		const interval& operator=(const double);

		const interval& operator+=(const interval&);
		const interval& operator-=(const interval&);
		const interval& operator*=(const double);
		const interval& operator/=(const double);
		const interval& operator*=(const interval&);
		const interval& operator/=(const interval&);

		void set_left (double l) { lft = l; }
		void set_right(double r) { rht = r; }

		double left() const{ return lft; }
		double right() const {return rht; }
	private:
		double lft;
		double rht;
	};

	/** interval non-member operators */

	const interval operator+(const interval& p, const interval& q);
	const interval operator-(const interval& p, const interval& q);
	const interval operator*(const interval& p, const double a);
	const interval operator*(const double a,	const interval& p);
	const interval operator*(const interval& p, const interval& q);
	const interval operator/(const interval& p, const interval& q);
	const interval operator/(const interval& p, const double a);

	std::ostream& operator<<(std::ostream& os, const interval& i);

}