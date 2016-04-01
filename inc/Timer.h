#ifndef TIMER_HPP
#define TIMER_HPP

#include <chrono>
#include <iostream>

/**
	Here we essentially write a wrapper for the chrono timers in C++ that are somewhat cumbersome
	to use in your programs.

	To use the timer class do the following:

	In our program we create a timer object and start "timing" by using the .start() function.
	After performing the operations we wish to time we finish "timing" using the .stop() function.
	We display the time in seconds to screen using the .display() function. This also returns
	the time duration as a double in case a user wishes to use or store this elsewhere. Simples.

	Note we have given three wrappers for each clock defined in chrono. However, on my system
	(Dell optiplex 7010, intel i7, Windows 7 professional) all three clocks are synonymous*. You
	will have to test to see if this is the case on your specific platform. From experience on my
	machine I would trust only the first three significant figures of the display; possibly the fourth. 
	Again this is something you should test for yourselves on your particular platform.

	*I might be misunderstanding the use of the clocks in chrono so feel free to correct my code.

	As a point about effective C++ code I have inherited from the chrono namespace, for example, the 
	timer class is a chrono::steady_clock class. But as I haven't written any of my own constructors 
	in the timer class default construction, copy construction, and copy-assignment is handled by 
	those functions in the base class, chrono::steady_clock in our example.
	Instead, we should use the "has a" relationship, i.e. the timer class should have a data member 
	that is a chrono::steady_clock class, rather than inheriting. 
*/
namespace phys{

	class timer : public std::chrono::steady_clock {
	public:
		void start(){  m_begin = std::chrono::steady_clock::now(); }
		void stop() {  m_end   = std::chrono::steady_clock::now(); }

		double display( bool supress = false ) const;

	private:
		using std::chrono::steady_clock::time_point;
		time_point m_begin;
		time_point m_end;
	};


	class hi_res_timer : public std::chrono::high_resolution_clock {

	public:
		void start(){  m_begin = std::chrono::high_resolution_clock::now(); }
		void stop() {  m_end   = std::chrono::high_resolution_clock::now(); }

		double display( bool supress = false ) const;

	private:
		using std::chrono::high_resolution_clock::time_point;
		time_point m_begin;
		time_point m_end;
	};


	class sys_timer : public std::chrono::system_clock {
	public:
		void start(){  m_begin = std::chrono::system_clock::now(); }
		void stop() {  m_end   = std::chrono::system_clock::now(); }

		double display( bool supress = false ) const;
	private:
		using std::chrono::system_clock::time_point;
		time_point m_begin;
		time_point m_end;
	};
}
#endif
