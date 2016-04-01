#include <adv/interval.hpp>

namespace phys{

	interval::interval() : lft(0.0), rht(0.0) {}
	interval::interval(double s, double e) : lft(s), rht(s) {}

	/** interval member operators */

	const interval& interval::operator=(const double a)
	{
		this->lft = 0;
		this->rht = a;
		return *this;
	}

	const interval& interval::operator+=(const interval& i)
	{
		this->lft += i.lft;
		this->rht += i.rht;
		return *this;
	}

	const interval& interval::operator-=(const interval& i)
	{
		this->lft -= i.rht;
		this->rht -= i.lft;
		return *this;
	}

	const interval& interval::operator*=(const double a)
	{
		this->lft *= a;
		this->rht *= a;
		return *this;
	}

	const interval& interval::operator/=(const double a)
	{
		this->lft /= a;
		this->rht /= a;
		return *this;
	}

	const interval& interval::operator*=(const interval& i)
	{
		using std::min;
		using std::max;
		this->lft = min(lft*i.lft, min(lft*i.rht, min(rht*i.lft, rht*i.rht)));
		this->rht = max(lft*i.lft, max(lft*i.rht, max(rht*i.lft, rht*i.rht)));
		return *this;
	}

	const interval& interval::operator/=(const interval& i)
	{
		using std::min;
		using std::max;
		//division by interval containing zero is not defined 
		if (i.lft * i.rht <= 0){
			std::cout << "Warning: dividing interval contains zero -> returning LHS.\n";
			return *this;
		}
		this->lft = min(lft / i.lft, min(lft / i.rht, min(rht / i.lft, rht / i.rht)));
		this->rht = max(lft / i.lft, max(lft / i.rht, max(rht / i.lft, rht / i.rht)));
		return *this;
	}

	/** interval non-member operators */

	const interval operator+(const interval& p, const interval& q)
	{
		return interval(p) += q;
	}

	const interval operator-(const interval& p, const interval& q)
	{
		return interval(p) -= q;
	}

	const interval operator*(const interval& p, const double a)
	{
		return interval(p) *= a;
	}

	const interval operator*(const double a, const interval& p)
	{
		return interval(p) *= a;
	}

	const interval operator*(const interval& p, const interval& q)
	{
		return interval(p) *= q;
	}

	const interval operator/(const interval& p, const interval& q)
	{
		return interval(p) /= q;
	}

	const interval operator/(const interval& p, const double a)
	{
		return interval(p) /= a;
	}

	std::ostream& operator<<(std::ostream& os, const interval& i)
	{
		os << "[" << i.left() << ", " << i.right() << "]\n";
		return os;
	}
}