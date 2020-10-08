#include "Storage.h"

#include "PhysicalUnits.h"

namespace phys {

//Class ODE_Storage implementation ------------------------------------------------

//Default constructor - use set(initial_system) after calling the default constructor
ODEStorage::ODEStorage() :
		m_independent(), m_dependent(), m_size_of_y(0), m_odeOrder(0), m_xName(
				"x"), m_yName("y"), m_dyName("dy") {
}

//Constructor using an (ODE) state - note the initial_system is not stored on construction
ODEStorage::ODEStorage(const phys::state &initial_system,
		const std::string &x_title, const std::string &y_title,
		const std::string &dy_title) :
		m_independent(), m_dependent(), m_xName(x_title), m_yName(y_title), m_dyName(
				dy_title) {
	set(initial_system); //sets m_size_of_y and m_odeOrder
}

void ODEStorage::set(const phys::state &system) {
	clear(); //ensures vectors empty before storing any new data.
	m_size_of_y = system.y.size(); // equals dims times ode order
	m_odeOrder = (system.dy.empty()) ? 1 : 2;
}

void ODEStorage::clear() {
	m_independent.clear();
	m_dependent.clear();
}

void ODEStorage::store(const phys::state &system) {
	m_independent.push_back(system.x);

	for (size_t i = 0; i < m_size_of_y; ++i) {
		m_dependent.push_back(system.y[i]);
	}
}

/* works even if vectors are (initially) empty */
void ODEStorage::append(const ODEStorage &add) {
	if (add.m_size_of_y != this->m_size_of_y && this->m_size_of_y > 0) {
		throw(std::runtime_error(
				"Additional data does not match the existing data - please check your code"));
	}
	m_independent.insert(m_independent.end(), add.m_independent.begin(),
			add.m_independent.end());
	m_dependent.insert(m_dependent.end(), add.m_dependent.begin(),
			add.m_dependent.end());
}

void ODEStorage::wrapTo2Pi(size_t dim) {
	/*	error check the choice of dimension wanted */
	if (dim * m_odeOrder > m_size_of_y) {
		__error_does_not_exist("Dependent Variable");
	}

	if (dim == 0) { //perform wrapping on all dimensions
		size_t numOfDims = m_size_of_y / m_odeOrder;
		for (size_t i = 0; i < m_dependent.size(); i += m_size_of_y) {
			for (size_t j = 0; j < numOfDims; ++j) {
				double temp = m_dependent[i + j];
				while (temp < -PI)
					temp += 2. * PI;
				while (temp >= +PI)
					temp -= 2. * PI;
				m_dependent[i + j] = temp;
			}
		}
	} else { //perform wrapping on requested dimension
		size_t base = dim - 1;
		for (size_t i = base; i < m_dependent.size(); i += m_size_of_y) {
			double temp = m_dependent[i];
			while (temp < -PI)
				temp += 2. * PI;
			while (temp >= +PI)
				temp -= 2. * PI;
			m_dependent[i] = temp;
		}
	}

}

void ODEStorage::write(const std::string &filename, bool write_headers) const {
	if (m_independent.empty()) {
		std::string errmsg = __FUNCTION__;
		errmsg +=
				" -- User Error: attempting to write data before any data stored.\n";
		throw(std::runtime_error(errmsg));
	}

	std::ofstream output_file(filename);

	if (output_file.is_open() == false) {
		std::cout << "Warning: " << filename
				<< " could not be opened for writing" << std::endl;
		return;
	}

	if (write_headers) {
		size_t loop_max = m_size_of_y / m_odeOrder;
		output_file << "#" << m_xName;
		for (size_t j = 0; j < loop_max; ++j) {
			output_file << "\t" << m_yName << j + 1;
			if (m_odeOrder == 2)
				output_file << "\t" << m_dyName << j + 1;
		}
		output_file << "\n";
	}

	for (size_t i = 0; i < m_independent.size(); ++i) {
		output_file << m_independent[i];
		for (size_t j = 0; j < m_size_of_y; ++j) {
			output_file << "\t" << m_dependent[m_size_of_y * i + j];
		}
		output_file << "\n";
	}

	output_file.close();
}

stdVec_d ODEStorage::get_independent() const {
	return m_independent;
}

stdVec_d ODEStorage::get_dependent(size_t dim, bool first_deriv) const {
	if (dim == 0) {
		//std::cout << "\nDimension zero requested; assuming you want dimension one\n";
		dim = 1;
	}

	/*	error check the choice of dimension wanted */
	if (dim * m_odeOrder > m_size_of_y) {
		if (first_deriv)
			__error_does_not_exist("First Derivative");
		else
			__error_does_not_exist("Dependent Variable");
	}

	/*	quick return for an ode of order 1 with 1 dimension*/
	if (m_size_of_y == 1)
		return m_dependent;

	size_t base = dim - 1;

	/*	If first derivative wanted offset base by num of dimensions
	 (see interleaving description of the ODE_Storage)	*/
	if (first_deriv)
		base += m_size_of_y / m_odeOrder;

	stdVec_d retval;
	for (size_t i = base; i < m_dependent.size(); i += m_size_of_y)
		retval.push_back(m_dependent[i]);

	return retval;
}

stdVec_d ODEStorage::get_first_deriv(size_t dim) const {
	/*	error check that we have an ode of order 2 */
	if (m_odeOrder != 2)
		__error_first_deriv();

	stdVec_d retval = get_dependent(dim, true);

	return retval;
}

void ODEStorage::__error_does_not_exist(std::string which) const {
	std::string errmsg = which;
	errmsg += ": the dimension chosen does not exist. Please check your code.";
	throw(std::runtime_error(errmsg));
}
void ODEStorage::__error_first_deriv() const {
	std::string errmsg =
			"This solution is for an ode of order 1, i.e. it does not have a\n";
	errmsg += "first derivative value computed. Please check your code.";
	throw(std::runtime_error(errmsg));
}

}//namespace
