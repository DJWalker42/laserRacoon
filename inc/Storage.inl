/* Implementation of the storage class*/

namespace phys {

template<typename T>
Storage<T>::Storage(const std::string &x_name, const std::string &y_name) :
		m_independent(), m_dependent(), m_multi(), m_xTitle(x_name), m_yTitle(
				y_name), m_keyTitles() {
}

template<typename T>
Storage<T>::Storage(const std::string &x_name,
		const std::vector<std::string> &key_names) :
		m_independent(), m_dependent(), m_multi(key_names.size()), m_xTitle(
				x_name), m_yTitle(), m_keyTitles(key_names) {
}

template<typename T>
void Storage<T>::clear() {
	m_independent.clear();
	m_dependent.clear();
	m_multi.clear();
}

template<typename T> inline
void Storage<T>::store(const T &x_val) {
	m_independent.push_back(x_val);
}

template<typename T> inline
void Storage<T>::store(const T &x_val, const T &y_val) {
	m_independent.push_back(x_val);
	m_dependent.push_back(y_val);
}

template<typename T> inline
void Storage<T>::store(const T &x_val, const T &y0_val, const T &y1_val) {
#if _DEBUG
			assert(m_multi.size() == 2);
#endif
	m_independent.push_back(x_val);
	m_multi[0].push_back(y0_val);
	m_multi[1].push_back(y1_val);
}

template<typename T> inline
void Storage<T>::store(const T &x_val, const std::vector<T> &y_vals) {
#if _DEBUG
			assert(m_multi.size() == y_vals.size());
#endif
	m_independent.push_back(x_val);
	for (size_t i = 0; i < y_vals.size(); ++i)
		m_multi[i].push_back(y_vals[i]);
}

template<typename T> inline
void Storage<T>::copy(const std::vector<T> &x_vals,
		const std::vector<T> &y_vals) {
	m_independent = x_vals;
	m_dependent = y_vals;
}

template<typename T> inline
void Storage<T>::copy(const std::vector<T> &x_vals,
		const std::vector<std::vector<T>> &y_vals) {
	m_dependent.clear();
	m_independent = x_vals;
	m_multi = y_vals;
}

template<typename T> inline
void Storage<T>::copy(const std::vector<T> &y_vals) {
	if (m_dependent.empty() == false) {
		m_multi.push_back(m_dependent);
		m_dependent.clear();
	}
	m_multi.push_back(y_vals);
}

template<typename T>
void Storage<T>::clear_names() {
	m_xTitle = "";
	m_yTitle = "";
	m_keyTitles.clear();
}

template<class T>
void Storage<T>::append(const Storage<T> &add) {
	m_independent.append(add.m_independent);
	if (!m_dependent.empty())
		m_dependent.append(add.m_dependent);
	if (!m_multi.empty())
		for (size_t i = 0; i < m_multi.size(); ++i)
			m_multi[i].append(add.get_dependent(i));
}

template<class T>
bool Storage<T>::read(const std::string &filename) {
	std::vector < T > data;
	size_t rows = 0, cols = 0;
	import_array(filename, data, rows, cols);
	if (data.empty())
		return false;
	if (cols == 1) {
		std::cout << "No m_dependent variable\n";
		return false;
	}
	if (cols > 2) // we have multiple m_dependent data
			{
		m_multi = std::vector < std::vector < T >> (cols - 1);
		for (size_t n = 0; n < rows; ++n) {
			m_independent.push_back(data[n * cols]);
			for (size_t i = 1; i < cols; ++i)
				m_multi[i - 1].push_back(data[n * cols + i]);
		}
	} else //cols == 2
	{
		for (size_t n = 0; n < rows; ++n) {
			m_independent.push_back(data[2 * n]);
			m_dependent.push_back(data[2 * n + 1]);
		}
	}
	return true;
}

template<class T>
void Storage<T>::write(const std::string &filename, bool headers) {
	if (m_multi.empty()) {
		if (write_single(filename, headers))
			return;
	} else {
		if (write_multi(filename, headers))
			return;
	}
	//If here then there has been an error opening a file for writing
	// -- display warning but don't throw error.
	std::cout << "Warning: " << filename << " could not be opened for writing"
			<< std::endl;
	return;
}

template<class T>
bool Storage<T>::write_single(const std::string &fn, bool hd) const {
	assert(m_independent.size() == m_dependent.size());
	std::ofstream output(fn);
	if (output.is_open() == false)
		return false;
	if (hd)
		output << "#" << m_xTitle << "\t" << m_yTitle << "\n";
	size_t n = m_independent.size();
	for (size_t i = 0; i < n; ++i)
		output << m_independent[i] << "\t" << m_dependent[i] << "\n";
	output.close();
	return true;
}

template<class T>
bool Storage<T>::write_multi(const std::string &fn, bool hd) {
	for (size_t i = 0; i < m_multi.size(); ++i)
		assert(m_independent.size() == m_multi[i].size());
	std::ofstream output(fn);
	if (output.is_open() == false)
		return false;

	if (hd) {
		if (m_keyTitles.empty()) {
			throw std::runtime_error("No key names given\n");
		}
		output << "#" << m_xTitle;
		for (size_t i = 0; i < m_keyTitles.size(); ++i)
			output << "\t" << m_keyTitles[i];
		output << "\n";
	}
	size_t n = m_independent.size(), m = m_multi.size();
	for (size_t j = 0; j < n; ++j) {
		output << m_independent[j];
		for (size_t i = 0; i < m; ++i)
			output << "\t" << m_multi[i][j];
		output << "\n";
	}
	output.close();
	return true;
}

} //namespace

