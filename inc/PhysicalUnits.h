#ifndef PHYSICAL_UNITS_HPP
#define PHYSICAL_UNITS_HPP

/*	
	Physical constants in SI units.
	Each unit is expressed as its significand only i.e. with a zero exponent. 
	The actual exponents are given in the comments.
*/

namespace phys {
	namespace constants {
		/* ----	Maths constants ---- */
		extern const double PI;			//!< ratio of a circle's circumference to its diameter (16 decimal places)
		extern const double s_p_m;		//!< seconds per minute (exp == e1)
		extern const double m_p_h;		//!< minutes per hour (exp == e1)
		extern const double h_p_d;		//!< hours per day (exp == e1)
		extern const double d_p_y;		//!< days per year (exp == e2) 

		/* ---- Cosmological Units ---- */
		extern const double UG;				//!< Universal Gravitational constant (exponent == e-11) m*m*m/kg/s/s
		extern const double solar_Mass;		//!< Mass of the sun (exponent == e30) kg
		extern const double yrs;			//!< yrs in seconds (3.15569)(exponent == e7) s								
		extern const double AU;				//!< astronomical units in m (exponent == e11) m
		extern const double sidereal_month; //!< one sidereal month in seconds(2.36059) (exponent == 6) s
		extern const double lunar_dist;		//!< distance earth to moon in metres (exponent == 8) m
		extern const double earth_mass;		//!< mass of the Earth in kg (exponent == 24) kg
		extern const double earth_radius;	//!< radius of the earth in metres (exponent == 6) [assuming a spherical earth] m
		extern const double moon_mass;		//!< the mass of the moon in kg (exponent == 22) kg
		extern const double moon_radius;	//!< radius of the moon in units of earth-moon distance (exponent == 6) m

		extern const double GS;	//!< gravitational constant in units of solar mass, yrs, and astronomical units.(39.4861) (no exponent).
		extern const double GE;	//!< gravitational constant in units of earth mass, sidereal month, and lunar distance(39.1021) (no exponent).

		/* ---- Quantum constants ---- */ 
		extern const double planck;				//!< Planck's constant (exponent == -34) m*m*kg/s OR J.s
		extern const double hbar;				//!< h/2PI (1.05457173)(exponent == -34) m*m*kg/s
		extern const double electron_mass;		//!< (rest) mass of an electron (exponent == -31) kg
		extern const double electron_charge;	//!< fundamental charge of an electron (exponent == -19) C
		extern const double planck_ev;			//!< Planck's constant in electron-volts (4.13566750) (exponent == -15) eV.s
		extern const double hbar_ev;			//!< hbar in electron-volts (6.58211926) (exponent == -16) eV.s
	}
}
#endif
