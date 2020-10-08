#include "Storage.h"

#include "PhysicalUnits.h"

namespace phys {

//Class ODE_Storage implementation ------------------------------------------------

//Default constructor - use set(initial_system) after calling the default constructor
ODEStorage::ODEStorage() :
		m_independent(), m_dependent(), m_size_of_y(0), m_odeOrder(0),
		m_xName("x"), m_yName("y"), m_dyName("dy") {
}

//Constructor using an (ODE) state - note the initial_system is not stored on construction
ODEStorage::ODEStorage(const phys::state &initial_system,
		const std::string &x_title, const std::string &y_title,
		const std::string &dy_title) :
		m_independent(), m_dependent(), m_xName(x_title), m_yName(y_title), m_dyName(
				dy_title) {
	set(initial_system); //sets m_size_of_y and m_odeOrder
}

ODEStorage::ODEStorage(const std::string& ode_store_file, int order, int dimensions)
{
	read(ode_store_file, order, dimensions);
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

	output_file << *this;

	output_file.close();
}

void ODEStorage::read(const std::string & filename, int ode_order, int num_dims) {
	if (ode_order > 2 || ode_order < 0) {
		throw std::runtime_error("1st or 2nd ode order only");
	}
	if (num_dims <= 0) {
		throw std::runtime_error("number of dimensions strictly positive required");
	}

	std::ifstream input_file(filename);
	if (input_file.is_open() == false) {
		throw std::runtime_error(std::string("Unable to open ") + filename);
	}

	m_odeOrder = ode_order;
	m_size_of_y = ode_order * num_dims;

	if(input_file.get() == '#') {
		//extract variable titles from the header
		std::string line;
		std::getline(input_file, line);
		std::stringstream line_stream(line);
		line_stream >> m_xName;
		line_stream >> m_yName;
		//remove the appended dimension number from the y (dy) name(s)
		m_yName.pop_back();
		if (ode_order == 2) {
			line_stream >> m_dyName;
			m_dyName.pop_back();
		}
	}

	//reset file to beginning
	input_file.seekg(0);

	input_file >> *this;
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

//make the ODE_Storage class serialise / deserialise using streams
std::ostream& operator<<(std::ostream &os, const ODEStorage &ode_s) {
	for (size_t i = 0; i < ode_s.m_independent.size(); ++i) {
		os << ode_s.m_independent[i];
		for (size_t j = 0; j < ode_s.m_size_of_y; ++j) {
			os<< "\t" << ode_s.m_dependent[ode_s.m_size_of_y * i + j];
		}
		os << "\n";
	}

	return os;
}
std::istream& operator>>(std::istream &is, ODEStorage &ode_s) {

	//clear any data
	ode_s.clear();

	//ignore the header line if present
	if (is.get() == '#') is.ignore(1000, '\n');

	/*
	 * Rows look like:
	 *  order 1: x y0 y1 y2 ....
	 *  order 2: x y0 dy0 y1 dy1 ...
	 *
	 *  numbers refer to dimensions in the system,
	 */
	std::string line;
	while (std::getline(is, line)){
		std::stringstream line_stream(line);
		int count{0};
		double value;
		while (line_stream >> value) {
			if (count++) {
				ode_s.m_dependent.push_back(value);
			} else {
				//zero case
				ode_s.m_independent.push_back(value);
			}
		}
	}

	//you will have to know the order of the ODE to which the data refers

	return is;
}




}//namespace
