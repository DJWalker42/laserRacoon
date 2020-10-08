#ifndef STORAGE_HPP
#define STORAGE_HPP

#define _CRT_SECURE_NO_WARNINGS

#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cassert>
#include <algorithm>

#include "State.h"

namespace phys {
/**	Class used to store the numerical solution of ODEs. As such all data stored will be of double type.
 There is no read function as it is assumed the data will come directly from the ODE solvers.
 If you wish to read in data use the general Storage class instead - but be aware that the
 general class has no idea about the interleaving strategy used by ODE_storage - the user
 will have to ensure they are getting the data they want if using the general storage.
 */
class ODEStorage {
public:
	/*	Default constructor - use set(initial_system) after calling the default constructor*/
	ODEStorage();
	/*	Constructor for Storage class using a differential state. Use the initial state of the system to
	 initialise a storage object; the data from the initial system used is NOT stored. The optional
	 string arguments provide headers for any output of the data.
	 */
	ODEStorage(const phys::state &initial_system,
			const std::string &x_title = "x", const std::string &y_title = "y",
			const std::string &dy_title = "dy");

	ODEStorage(const std::string& ode_store_file, int order, int dimensions);

public: //interface functions
	/*	Function to clear the storage vectors
	 */
	void clear();
	/* Function set/reset a storage object to a new differential system
	 */
	void set(const phys::state &system);
	/*	Function to store the system data.
	 If default constructor used then use .set(initial_system) to prepare storage class for data.
	 The independent variable is stored to the vector independent_var.
	 The dependent variable(s) are stored to the vector dependent_vars in the following manner:

	 For a 1st order system the dependent vector will contain the interleaved variables
	 from each dimension. For example, a 1st order system with 3 dimensions (x,y,z) is stored as
	 dependent_vars[0] = x0, dependent_vars[1] = y0, dependent_vars[2] = z0, dependent_vars[3] = x1,
	 and so on, where the number after the variable refers to the integration step index.

	 For a 2nd order system the dependent vector will contain the variables (e.g y) up to the dimensions,
	 then the first derivative of those variables (e.g. dy), up to the dimensions.
	 For example, a 2nd order system with 2 dimensions (x,y) is stored as dependent_vars[0] = x0,
	 dependent_vars[1] = y0, dependent_vars[2] = dx0 , dependent_vars[3] = dy0, dependent_vars[4] = x1,
	 dependent_vars[5] = y1, and so on.
	 */
	void store(const phys::state &system);

	/*	Appends data from another ODE_Storage object to *this Storage object
	 Warning -- if the data being appended is not of the same dimensions & order of the existing data
	 will throw an exception */
	void append(const ODEStorage &additional);

	/*	Function to wrap the dependent data post solution to either -pi to +pi. The user may
	 select a particular dimension to wrap; the default is to wrap all dimensions.

	 Note: This function works by adding or subtracting 2PI until the value falls within the specified range.
	 Therefore the error in the wrapped values is equal to +/- N times the machine epsilon for doubles where N is
	 the number of additions/subtractions required. In extreme cases this may lead to values not being wrapped correctly.
	 */
	void wrapTo2Pi(size_t dimension = 0);

	/*	Function to set names of variables -- option to write these as headers of save file.
	 They are also used by the graphing functions as titles for axes and keys.
	 */
	void set_x_name(const std::string &name) {
		m_xName = name;
	}
	void set_y_name(const std::string &name) {
		m_yName = name;
	}
	void set_dy_name(const std::string &name) {
		m_dyName = name;
	}

	/*	Functions to get data names	*/
	std::string get_x_name() const {
		return m_xName;
	}
	std::string get_y_name() const {
		return m_yName;
	}
	std::string get_dy_name() const {
		return m_dyName;
	}

	/*Functions to get size information. */
	size_t num_grid_points() const {
		return m_independent.size();
	}
	size_t num_of_vars() const {
		return m_dependent.size();
	} // should equal order * dims * num_grid_points

	/*	Function to write the stored system to <filename>. Can handle both first and second ordered systems
	 with multiple dimensions.*/
	void write(const std::string &filename, bool write_headers = false) const;

	/*
	 * Function to read ODEStorage data from the provided file, file created by the 'write'
	 * member function
	 */
	void read(const std::string & filename, int ode_order, int num_dims);

	/*	Returns the stepped independent variable (i.e. start -> end in steps of h)
	 Use after applying a full_solve() to get access to the data.*/
	stdVec_d get_independent() const;

	/*	Returns the solved dependent variable of the dimension of choice, dim.
	 Default dim = 1. Use after applying a full_solve() to get access to the data.*/
	stdVec_d get_dependent(size_t dim = 1, bool first_deriv = false) const;

	/*	Returns the solved first derivative of the dependent variable of the dimension of choice, dim.
	 Default dim = 1. Use after applying a full_solve() to get access to the data.*/
	stdVec_d get_first_deriv(size_t dim = 1) const;

	//make the ODE_Storage class serialise / deserialise using streams
	friend std::ostream& operator<<(std::ostream &os, const ODEStorage &ode_s);
	friend std::istream& operator>>(std::istream &is, ODEStorage &ode_s);

private: //functions
	void __error_does_not_exist(std::string which) const;
	void __error_first_deriv() const;
private: // internal representation
	stdVec_d m_independent;
	stdVec_d m_dependent;
	size_t m_size_of_y;
	int m_odeOrder;
	std::string m_xName;
	std::string m_yName;
	std::string m_dyName;
};

/* -------------------------------------------------------------------------------------------*/

/** --- Storage for any other type of data that has an independent and dependent values --- */
template<class T>
class Storage {
public: //interface
	/* Default constructor */
	Storage(const std::string &xName = "x", const std::string &yName = "y");
	/* Constructor to set names for multiple dependent data */
	Storage(const std::string &xName, const std::vector<std::string> &yNames);

	//clears all data
	void clear();
	//stores a single value of the independent variable
	inline void store(const T &x_val);
	//stores a single value of both the independent and dependent variables.
	inline void store(const T &x_val, const T &y_val);
	//stores a single value of the independent and two dependent variables. m_multi must have size 2.
	inline void store(const T &x_val, const T &y0_val, const T &y1_val);
	//stores a single value of the independent and multiple dependent values. y_vals size must match m_multi size.
	inline void store(const T &x_val, const std::vector<T> &y_vals);

	//stores vectors of values of both the independent and dependent variables.
	inline void copy(const std::vector<T> &x_vals,
			const std::vector<T> &y_vals);
	//stores a vector of the independent value, and multiple vectors of dependent variables. Clears any dependent values.
	inline void copy(const std::vector<T> &x_vals,
			const std::vector<std::vector<T>> &y_vals);
	//stores a single vector of dependent variables -- if the dependent vector already has
	//data, assumes you want to store multiple dependent variables; does not overwrite data.
	inline void copy(const std::vector<T> &y_vals);

	void set_x_name(const std::string &name) {
		m_xTitle = name;
	}
	void set_key_names(const std::vector<std::string> &names) {
		m_keyTitles = names;
	}
	void set_y_name(const std::string &name) {
		m_yTitle = name;
	}
	void add_to_key(const std::string &name) {
		m_keyTitles.push_back(name);
	}

	void clear_names();

	//This function appends data from another storage object ('add') to *this storage object
	//Assumes additional data multi member has same dimensions as the existing multi member
	//- if the additional has less an exception will be thrown, if it has more the data in the
	//extra vectors will be ignored.
	void append(const Storage<T> &add);

	/** Read in data from a text file. Format of the file should be
	 x1 a1 b1 c1 ....
	 x2 a2 b2 c2 ....
	 and so on
	 where x is the independent variable and a,b,c,... are the dependent variables measured
	 or computed at x.
	 */
	bool read(const std::string &filename);

	//writes stored data to a tab delimited text file (x	y1	[z1 ...] )
	//can handle both single and multiple data; if headers wanted set second argument to true.
	void write(const std::string &filename, bool headers = false);

	//Query data
	std::vector<T> get_independent() const {
		return m_independent;
	}

	std::vector<T> get_dependent() const {
		return m_dependent.empty() ? get_dependent(0) : m_dependent;
	}

	std::vector<std::vector<T>> get_multi() const {
		return m_multi;
	}
	std::vector<T> get_dependent(size_t dim) const {
		if (m_multi.empty())
			throw(std::runtime_error("No multi-data stored"));
		if (dim > m_multi.size() - 1)
			throw(std::runtime_error("Dimension too large"));
		return m_multi[dim];
	}
	//Query names
	std::string get_x_name() const {
		return m_xTitle;
	}
	std::string get_y_name() const {
		return m_yTitle;
	}
	std::string get_name(size_t index) const {
		assert(index < m_keyTitles.size());
		return m_keyTitles[index];
	}
	std::vector<std::string> get_key_names() const {
		return m_keyTitles;
	}
private: //functions
	bool write_single(const std::string &fn, bool hd) const;
	bool write_multi(const std::string &fn, bool hd);
private: // representation
	std::vector<T> m_independent;			//!< e.g. x
	std::vector<T> m_dependent;				//!< e.g. y (= f(x))
	std::vector<std::vector<T>> m_multi;//!< multiple dependent variables for one independent variable
	std::string m_xTitle;		//!< For graphical purposes; name of the x-data.
	std::string m_yTitle;		//!< For graphical purposes; name of the y-data.
	std::vector<std::string> m_keyTitles;//!< For graphical purposes; name of the multiple y-data
};

} //namespace

#include "Storage.inl"

#endif //header guard
