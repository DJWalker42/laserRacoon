#include <iostream>
#include <cstdint>
#include <numeric>

#include "Timer.h"
#include "Storage.h"
#include "Visualise.h"

/*
 *  We need to assume a format for the Matrix classes;
 *  we'll assume row major format in memory layout
 *
 * Matrix_t1
 *  By using a data pointer with a single level of indirection
 *  to represent a two-dimensional matrix we have to compute a
 *  linear index for a given coordinate (i,j); this adds extra
 *  floating point operations into our code, but which are
 *  likely removed by automatic compiler optimisations, such as
 *  loop unrolling (using -O1 or higher)
 */
class Matrix_t1 {
private:
	double * _data;
	uint32_t _rows;
	uint32_t _cols;
public:
	Matrix_t1(uint32_t rows, uint32_t cols, double val = 0.) : _rows(rows), _cols(cols)
	{
		assert(rows != 0 && cols != 0);

		_data = new double [_rows * _cols];
		for (int i = 0; i < _rows * _cols; ++i) {
			_data[i] = val;
		}

	}

	Matrix_t1(const Matrix_t1& other) : Matrix_t1(other._rows, other._cols)
	{
		//Do we need to assert that rows and columns are strictly positive here?

		int N = _rows * _cols;
		for (int i = 0; i < N; ++i) {
			_data[i] = other._data[i];
		}
	}

	~Matrix_t1()
	{
		delete[] _data;
	}

	const Matrix_t1& operator=(const Matrix_t1& other) {
		//check for self-assignment
		if (this == &other) return *this;

		//Do we need to assert that rows and columns are strictly positive here?
		//Note that to use an assignment operator, *this object must be initialised
		//(constructed) before hand

		 _rows = other._rows;
		 _cols = other._cols;

		int N = _rows * _cols;
		for (int i = 0; i < N; ++i) {
			_data[i] = other._data[i];
		}

		return *this;
	}

	/*
	 *  By defining a copy constructor and an copy assignment operator we suppress
	 *  default (compiler defined) move versions of those functions.
	 */

public:
	//element read access
	double operator()(uint32_t i, uint32_t j) const {
		return _data[j + _cols * i];
	}

	//element write access
	double& operator()(uint32_t i, uint32_t j) {
		return _data[j + _cols * i];
	}

	Matrix_t1 operator*(const Matrix_t1& rhs) {

		assert(this->_cols == rhs._rows);

		Matrix_t1 result(this->_rows, rhs._cols);

		//naive matrix multiplication
		for (int i = 0; i < this->_rows; ++i) {
			for (int j = 0; j < rhs._cols; ++j) {
				for (int k = 0; k < this->_cols; ++ k) {
					result._data[i + j * this->_rows] +=
							this->_data[k + i * this->_cols] * rhs._data[j + k * rhs._cols];
				}
			}
		}

		return result;
	}

	//prints as [row]\n[row]\n ... up to number of rows.
	friend std::ostream & operator<<(std::ostream& os, const Matrix_t1 & mat) {
		for (int i = 0; i < mat._rows; ++i) {
			for (int j = 0; j < mat._cols; ++j) {
				os << mat._data[j + i * mat._cols] << "\t";
			}
			os << "\n";
		}
		os << "\n\n";
		return os;
	}
}; //end Matrix_t1

/*
 *  Matrix_t2 class uses a vector of vectors to represent a matrix.
 *  Note that the "head" of a vector is stored in the stack but its
 *  data is stored on the heap. You should also find out the difference
 *  between the size of a vector and its capacity, and how those
 *  mechanisms could affect the code's performance.
 */
class Matrix_t2 {
private:
	std::vector<std::vector<double>> _data;
public:
	Matrix_t2(uint32_t rows, uint32_t cols, double val = 0.) :
		_data(rows, std::vector<double>(cols, val))
	{
		assert(rows != 0 && cols != 0);
	}

	/*
	 *  The default copy constructor and assignment operator do what
	 *  we want due to the use of vectors as the internal representation.
	 */

public:
	//element read access
	double operator()(uint32_t i, uint32_t j) const {
		return _data[i][j];
	}

	//element write access
	double& operator()(uint32_t i, uint32_t j) {
		return _data[i][j];
	}

	//row read access
	std::vector<double> operator()(uint32_t i) const {
		return _data[i];
	}

	//row write access (why might we not want to give write access to a row, which is a std::vector type?)
	std::vector<double>& operator()(uint32_t i) {
		return _data[i];
	}

	uint32_t rows() const noexcept { return _data.size(); }
	uint32_t cols() const noexcept { return _data[0].size(); }


public:

	Matrix_t2 operator*(const Matrix_t2& rhs) {
		uint32_t Arows = _data.size();
		uint32_t Acols = _data[0].size();
		uint32_t Brows = rhs._data.size();
		uint32_t Bcols = rhs._data[0].size();

		assert(Acols == Brows);

		Matrix_t2 result(Arows, Bcols);

		//naive matrix multiplication
		for (int i = 0; i < Arows; ++i) {
			for (int j = 0; j < Bcols; ++j) {
				for (int k = 0; k < Acols; ++ k) {
					result._data[i][j] += _data[i][k] * rhs._data[k][j];
				}
			}
		}
		return result;
	}

	friend std::ostream& operator<<(std::ostream& os, const Matrix_t2& mat) {
		for (int i = 0; i < mat._data.size(); ++i) {
			for (int j = 0; j < mat._data[0].size(); ++j) {
				os << mat._data[i][j] << "\t";
			}
			os << "\n";
		}
		os << "\n\n";
		return os;
	}
}; //end Matrix_t2

/*
 *  Matrix_t3 class uses two levels of indirection to model a
 *  two-dimensional matrix.
 *  Note that the memory layout of the internal representation
 *  of this matrix will not be in a contiguous block, unlike making
 *  a static array using the syntax:
 *
 *  double _data [rows][cols];
 *
 *  but note that 'rows' and 'cols' need to be compile time
 *  constants or constant expressions in this case.
 */
class Matrix_t3 {
private:
	double ** _data;
	uint32_t _rows;
	uint32_t _cols;
public:
	Matrix_t3(uint32_t rows, uint32_t cols, double val = 0.) :
		_rows(rows), _cols(cols) {
		assert(_rows != 0 && _cols != 0);

		_data = new double* [_rows];
		for (int i = 0; i < _rows; ++i) {
			_data[i] = new double[_cols];
			for (int j = 0; j <_cols; ++j) {
				_data[i][j] = val;
			}
		}
	}

	~Matrix_t3() {
		//delete nested memory first
		for (int i = 0; i < _rows; ++i) {
			delete [] _data[i];
		}
		delete [] _data;
	}

	Matrix_t3(const Matrix_t3& other) : Matrix_t3(other._rows, other._cols)
	{
		for (int i = 0; i < _rows; ++i) {
			for (int j = 0; j <_cols; ++j) {
				_data[i][j] = other._data[i][j];
			}
		}
	}

	const Matrix_t3& operator=(const Matrix_t3& other) {
		//check for self-assignment
		if (this == &other) return *this;

		 _rows = other._rows;
		 _cols = other._cols;

		for (int i = 0; i < _rows; ++i) {
			for (int j = 0; j <_cols; ++j) {
				_data[i][j] = other._data[i][j];
			}
		}

		return *this;
	}

	//element read access
	double operator()(uint32_t i, uint32_t j) const {
		return _data[i][j];
	}

	//element write access
	double& operator()(uint32_t i, uint32_t j) {
		return _data[i][j];
	}

	Matrix_t3 operator*(const Matrix_t3& rhs) {

		assert(this->_cols == rhs._rows);

		Matrix_t3 result(this->_rows, rhs._cols);

		//naive matrix multiplication
		for (int i = 0; i < this->_rows; ++i) {
			for (int j = 0; j < rhs._cols; ++j) {
				for (int k = 0; k < this->_cols; ++ k) {
					result._data[i][j] += _data[i][k] * rhs._data[k][j];
				}
			}
		}

		return result;
	}

	//prints as [row]\n[row]\n ... up to number of rows.
	friend std::ostream & operator<<(std::ostream& os, const Matrix_t3 & mat) {
		for (int i = 0; i < mat._rows; ++i) {
			for (int j = 0; j < mat._cols; ++j) {
				os << mat._data[i][j] << "\t";
			}
			os << "\n";
		}
		os << "\n\n";
		return os;
	}
}; //end Matrix_t3



int main () {

	const int max_repeats {10};
	const double processor_freq {2.2e9}; //my machine 2.2GHz processor (may need to change)

	phys::timer t;
	std::vector<std::string> names {"Model", "Matrix_t1", "Matrix_t2", "Matrix_t3"};
	phys::Storage<double> container("n", names);

	for (int n = 150; n <= 500; n += 1) {
		Matrix_t1 A(n, n, 1.);
		Matrix_t1 B(n, n, 1.);

		Matrix_t2 D(n, n, 1.);
		Matrix_t2 E(n, n, 1.);

		Matrix_t3 G(n, n, 1.);
		Matrix_t3 H(n, n, 1.);

		std::vector<double> times;

		//model: given A(m,n) B(n,p), A * B takes 2mnp - mp flops => square 2n^3 - n^2
		// the -n^2 is insignificant for n >> 1 but leave in for accuracy.
		int n2 = n * n;
		times.push_back((2. * n2 * n - n2)/processor_freq);

		t.start();
		for (int repeat = 0; repeat < max_repeats; ++repeat) {
			Matrix_t1 C  = A * B;
		}
		t.stop();

		times.push_back(t.get() / max_repeats);


		t.start();
		for (int repeat = 0; repeat < max_repeats; ++repeat) {
			Matrix_t2 F  = D * E;
		}
		t.stop();

		times.push_back(t.get() / max_repeats);

		t.start();
		for (int repeat = 0; repeat < max_repeats; ++repeat) {
			Matrix_t3 P  = G * H;
		}
		t.stop();

		times.push_back(t.get() / max_repeats);

		container.store(double(n), times);

		std::cout << "n = " << n << std::endl;
 	}

	container.write("./matrix_multiply.log", true);

	phys::Viewer viewer(1600, 900);
	viewer.plot(container);

	return 0;

}
