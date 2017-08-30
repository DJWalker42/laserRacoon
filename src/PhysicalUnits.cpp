#include "PhysicalUnits.h"

/*	Any "magic" numbers found in expressions, e.g. * 10.0 or * 1e-2, 
	are required to obtain the correct exponent calculation. 

	For example,

	GE is computed from the following:
	UG = 6.67384e-11 mmm/kg/s/s (universal gravitational constant)
	mass_unit = 5.97219e24 (kg) (earth's mass)
	time_unit = 2.36059e6 (s) (sidereal month in seconds)
	distance_unit = 3.844e8 (m) (lunar distance; earth to moon)

	GE = UG * mass_unit * time_unit * time_unit / distance_unit / distance_unit / distance_unit
	GE = (6.67384 * 5.97219 * 2.36059 * 2.36059 / 3.844 / 3.844 / 3.844) * pow(10, -11 + 24 + 6 + 6 - 8 - 8 - 8)
	GE = 3.9102 * pow(10, 1) = 3.9102 * 10

	hence the * 10 in the expression.
*/

namespace phys{
	namespace constants{
		/* ----	Maths constants ---- */
		const double PI{3.1415926535897932e0};
		const double s_p_m{6.e0};
		const double m_p_h{6.e0};
		const double h_p_d{2.4e0};
		const double d_p_y{3.65242e0};

		/* ---- Cosmological Units ---- */
		const double UG{6.67384e0};
		const double solar_Mass{1.9891e0};
		const double yrs{s_p_m * m_p_h * h_p_d * d_p_y * 1e-2};
		const double AU{1.495978707e0};
		const double sidereal_month{2.7321661 * s_p_m * m_p_h * h_p_d * 1e-2};
		const double lunar_dist{3.844e0};
		const double earth_mass{5.97219e0};
		const double earth_radius{6.371e0};
		const double moon_mass{7.34767309e0};
		const double moon_radius{1.7371e0};

		const double GS{UG * solar_Mass * yrs * yrs / AU / AU / AU};
		const double GE{UG * 10 * earth_mass * sidereal_month * sidereal_month / lunar_dist / lunar_dist / lunar_dist};

		/* ---- Quantum constants ---- */
		const double planck{6.62606957e0};
		const double hbar{planck / 2 / PI};
		const double electron_mass{9.10938291e0};
		const double electron_charge{1.60217657e0};
		const double planck_ev{planck / electron_charge};
		const double hbar_ev{10.0*hbar / electron_charge};

	}
}
