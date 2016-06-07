#include "PhysicalUnits.h"

namespace phys{
	namespace constants{
		/* ----	Maths constants ---- */
		extern const double PI{3.1415926535897932e0};	
		extern const double s_p_m{6.e0};				
		extern const double m_p_h{6.e0};				
		extern const double h_p_d{2.4e0};			
		extern const double d_p_y{3.65242e0};		 

		/* ---- Cosmological Units ---- */
		extern const double UG{6.67384e0};				
		extern const double solar_Mass{1.9891e0};			
		extern const double yrs{s_p_m * m_p_h * h_p_d * d_p_y * 1e-2};
		extern const double AU{1.495978707e0};			
		extern const double sidereal_month{2.7321661 * s_p_m * m_p_h * h_p_d * 1e-2};			
		extern const double lunar_dist{3.844e0};			
		extern const double earth_mass{5.97219e0};		
		extern const double earth_radius{6.371e0};		
		extern const double moon_mass{7.34767309e0};		
		extern const double moon_radius{1.7371e0};		

		extern const double GS{UG * solar_Mass * yrs * yrs / AU / AU / AU};	
		extern const double GE{UG * 10 * earth_mass * sidereal_month * sidereal_month / lunar_dist / lunar_dist / lunar_dist};	

		/* ---- Quantum constants ---- */
		extern const double planck{6.62606957e0};					
		extern const double hbar{planck / 2 / PI};					
		extern const double electron_mass{9.10938291e0};			
		extern const double electron_charge{1.60217657e0};		
		extern const double planck_ev{planck / electron_charge};	
		extern const double hbar_ev{10.0*hbar / electron_charge};	

	}
}
