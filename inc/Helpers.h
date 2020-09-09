#ifndef HELPERS_HPP
#define HELPERS_HPP

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

namespace phys{


/*toUType - function to return the underlying type of an enumerator
 * Typically used with scoped (class) enums but can be used with any enum structure
 */
template<typename E>
constexpr typename std::underlying_type<E>::type
toUType(E enumerator) noexcept{
  return static_cast<typename std::underlying_type<E>::type> (enumerator);
}



/*	Read in a single vector by reference from a string of numbers; use getline to get string.
	Returns size of the vector */
template<class T>
size_t import_vector(const std::string& line, std::vector<T>& vec)
{
	std::istringstream is(line);
	T n;
	while (is >> n)
		vec.push_back(n);
	return vec.size();
}

/*	Read in an array of numbers from text file.
	Appends values to the vector vec*/
template<class T>
void import_array(const std::string& filename, std::vector<T>&vec, size_t& rows, size_t& cols)
{
	std::ifstream file(filename);
	if(file.is_open())
	{
		std::string line;
		std::getline(file, line);
		++rows;
		cols = import_vector(line, vec);

		while(std::getline(file, line))
		{
			++rows;
			import_vector(line, vec);
		}
		file.close();
	}
	else
	{
		std::cout << "Error opening " << filename << "\n";
		std::cout << "Vector argument unmodified\n";
	}
}

template<class T>
void swap( T& a, T& b)
{
	T c = a;
	a = b;
	b = c;
}

} //namespace

#endif //header guard
