#include "Timer.h"

namespace phys {

void timer::display() const {
	using namespace std::chrono;
	if (m_end < m_begin) {
		std::cout << "Did you stop() before start()?\n";
		return;
	}
	duration<double> time_span = duration_cast<duration<double>>(m_end - m_begin);

	std::cout << "Execution time: " << time_span.count() << " seconds.\n";
}


} //namespace
