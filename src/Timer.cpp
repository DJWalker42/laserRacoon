#include "Timer.h"

namespace phys{

	typedef std::chrono::duration<double> duration_d;


	double timer::display( bool sup ) const
	{
		using namespace std::chrono;
		if(m_end < m_begin)
		{
			std::cout << "Did you stop() before start()?\n";
			return 0.;
		}
		duration_d time_span = duration_cast<duration_d>(m_end - m_begin);
		if(!sup)
			std::cout << "Execution time: " << time_span.count() << " seconds.\n";
		return time_span.count();
	}

	double hi_res_timer::display( bool sup ) const
	{
		using namespace std::chrono;
		if(m_end < m_begin)
		{
			std::cout << "Did you stop() before start()?\n";
			return 0.;
		}
		duration_d time_span = duration_cast<duration_d>(m_end - m_begin);
		if(!sup)
			std::cout << "Execution time: " << time_span.count() << " seconds.\n";
		return time_span.count();
	}

	double sys_timer::display( bool sup ) const
	{
		using namespace std::chrono;
		if(m_end < m_begin)
		{
			std::cout << "Did you stop() before start()?\n";
			return 0.;
		}
		duration_d time_span = duration_cast<duration_d>(m_end - m_begin);
		if(!sup)
			std::cout << "Execution time: " << time_span.count() << " seconds.\n";
		return time_span.count();
	}

}
