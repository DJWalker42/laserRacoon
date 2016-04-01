#ifndef PHYSICAL_UNITS_HPP
#define PHYSICAL_UNITS_HPP

/*	
	Physical constants in SI units.
	Each unit is expressed as its mantissa only. The actual exponents are given in the comments
*/

/* ----	Maths constants ---- */
const double PI = 3.1415926535897932e0;	//!< ratio of a circle's circumference to its diameter (16 decimal places)
const double s_p_m = 6.e0;				//!< seconds per minute (exp == e1)
const double m_p_h = 6.e0;				//!< minutes per hour (exp == e1)
const double h_p_d  = 2.4e0;			//!< hours per day (exp == e1)
const double d_p_y  = 3.65242e0;		//!< days per year (exp == e2) 

/* ---- Cosmological Units ---- */
const double UG = 6.67384e0;				//!< Universal Gravitational constant (exponent == e-11) m*m*m/kg/s/s
const double solar_Mass = 1.9891e0;			//!< Mass of the sun (exponent == e30) kg
const double yrs = s_p_m * m_p_h *			//!< yrs in seconds (3.15569)(exponent == e7) s
	h_p_d * d_p_y * 1e-2;					
const double AU = 1.495978707e0;			//!< astronomical units in m (exponent == e11) m
const double sidereal_month = 2.7321661 * 
	s_p_m * m_p_h * h_p_d * 1e-2;			//!< one sidereal month in seconds(2.36059) (exponent == 6) s
const double lunar_dist = 3.844e0;			//!< distance earth to moon in metres (exponent == 8) m
const double earth_mass = 5.97219e0;		//!< mass of the Earth in kg (exponent == 24) kg
const double earth_radius = 6.371e0;		//!< radius of the earth in metres (exponent == 6) [assuming a spherical earth] m
const double moon_mass = 7.34767309e0;		//!< the mass of the moon in kg (exponent == 22) kg
const double moon_radius = 1.7371e0;		//!< radius of the moon in units of earth-moon distance (exponent == 6) m

const double GS = UG * solar_Mass * yrs * yrs /AU/AU/AU;	//!< gravitational constant in units of solar mass, yrs, and astronimcal units.(39.4861) (no exponent).
const double GE = UG * 10 * earth_mass * sidereal_month * 
	sidereal_month / lunar_dist /lunar_dist/ lunar_dist;	//!< gravitational constant in units of earth mass, sidereal month, and lunar distance(39.1021) (no exponent).

/* ---- Quantum constants ---- */ 
const double planck = 6.62606957e0;					//!< Planck's constant (exponent == -34) m*m*kg/s OR J.s
const double hbar = planck/2/PI;					//!< h/2PI (1.05457173)(exponent == -34) m*m*kg/s
const double electron_mass = 9.10938291e0;			//!< (rest) mass of an electron (exponent == -31) kg
const double electron_charge = 1.60217657e0;		//!< fundamental charge of an electron (exponent == e-19) C
const double planck_ev = planck/electron_charge;	//!< plancks constant in electron-volts (4.13566750) (expoenent == -15) eV.s
const double hbar_ev = 10.0*hbar/electron_charge;	//!< hbar in electron-volts (6.58211926) (expoenent == -16) eV.s
#endif
