#ifndef TIMER_HPP
#define TIMER_HPP

#include <chrono>
#include <iostream>

/**
 Here we essentially write a wrapper for a chrono timers in C++ that are somewhat cumbersome
 to use in your programs.

 To use the timer class do the following:

 In our program we create a timer object and start "timing" by using the .start() function.
 After performing the operations we wish to time we finish "timing" using the .stop() function.
 We display the time in seconds to screen using the .display() function. If you just want the
 time in seconds use the get() member function

 As a point about effective C++ code I have inherited from the chrono namespace, for example, the
 timer class is a chrono::steady_clock class. But as I haven't written any of my own constructors
 in the timer class default construction, copy construction, and copy-assignment is handled by
 those functions in the base class, chrono::steady_clock in our example.
 Instead, we should use the "has a" relationship, i.e. the timer class should have a data member
 that is a chrono::steady_clock class, rather than inheriting.
 */
namespace phys {

class timer {
public:
	void start() {
		m_begin = std::chrono::steady_clock::now();
	}
	void stop() {
		m_end = std::chrono::steady_clock::now();
	}

	void display() const;

	double get() const noexcept {
		using namespace std::chrono;
		duration<double> time_span = duration_cast<duration<double>>(m_end - m_begin);
		return time_span.count();
	}

private:
	using time_point = std::chrono::steady_clock::time_point;
	time_point m_begin;
	time_point m_end;
};


} //namespace

#endif //header guard
