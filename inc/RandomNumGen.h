#ifndef DJW_RANDOMNUMGEN_HPP
#define DJW_RANDOMNUMGEN_HPP

#include <random>
#include <chrono>

#include "DynVector.h"

/*
 Here we note down all possible distributions for the random number
 generators available in the C++ library. To template for all
 distributions we note down their parameter(s) input types and
 output types. The user can use this list to select the specific
 constructor to use for a particular distribution.

 ----------------------------------------------------------------------------------
 input1	[input2]	output
 bernoulli_distribution				double	 NULL		bool
 binomial_distribution				int_t	 double		int_t   -- separate class
 cauchy_distribution					real_t	 real_t		real_t
 chi_squared_distribution			real_t	 NULL		real_t
 discrete_distribution				***	special			int_t
 exponential_distribution			real_t	 NULL		real_t
 extreme_value_distribution			real_t	 real_t		real_t
 fisher_f_distribution				real_t	 real_t		real_t
 gamma_distribution					real_t	 real_t		real_t
 geometric_distribution				double	 NULL		int_t
 lognormal_distribution				real_t	 real_t		real_t
 negative_binomial_distribution		int_t	 double		int_t  -- separate class
 normal_distribution					real_t	 real_t		real_t
 piecewise_constant_distribution		***	special			real_t
 piecewise_linear_distribution		***	special			real_t
 poisson_distribution				double	 NULL		int_t
 student_t_distribution				real_t	 NULL		real_t
 uniform_int_distribution			int_t	 int_t		int_t
 uniform_real_distribution			real_t	 real_t		real_t
 weibull_distribution				real_t	 real_t		real_t

 *** special
 size_t num_weights, double xmin, double xmax, fw (with xmax > xmin)
 where fw is a function pointer or object function that takes a double and returns a double
 fw is passed the values (xmin + (xmax- xmin)*(k+0.5)/nw), k = 0 -> nw-1,
 nw == num of weights (discrete)
 or num of subintervals (peicewise)


 ---------------------------------------------------------------------------------------------

 Examples of use: (std::minstd_rand is the default engine; time count is the default seed)

 phys::monte::RNG<> default_generator(0.0, 1.0)
 //uniform distribution on [0,1) returning doubles

 phys::monte::RNG<	double,
 std::minstd_rand,
 std::bernoulli_distribution,
 bool >	bernoulli_generator(0.5);
 //bernoulli distribution with a probability of 0.5.
 //Note that the bernoulli distribution is *NOT* templated; it *HAS* to return a bool.

 phys::monte::RNG<	std::minstd_rand,
 std::geometric_distribution<int>,
 int,
 double >	geomertric_generator(0.3)
 //geometric distribution with a 0.3 probability of success returing int type.

 phys::monte::RNG <	std::minstd_rand,
 std::discrete_distribution<int>,
 int >	discrete_generator(4, 0.0, 40.0, plus5)
 //discrete generator that will give the probabilities {0.1, 0.2, 0.3, 0.4}
 //plus5 is  function that adds 5 to a value and returns the result
 //values passed to the function are computed using (xmin + (xmax- xmin)*(k+0.5)/n)
 //note that for discrete, peicewise_constant, and peicewise_linear distributions
 //inputType is redundant.

 General template filler form:
 phys::monte::RNG<	input1_type,
 choice_of_engine,
 distribution<input1_type>,
 output_type
 >

 Note that both the binomial and negative binomial have been seperated due to ambiguities
 in the constructor and member function signatures when templated.

 phys::monte::RNG_binomial<> binomial_generator(3, 0.6)
 // binomial distribution with upper bound of 3 and probablity 0.6, returning int types.

 For more details go to:
 http://www.cplusplus.com/reference/random/

 TODO:	1. Change the system_clock to high_resolution_clock on the unsigned seeds
 2. Change functional style casts to static_casts


 */
namespace phys {
template<class inputType = double, class rngType = std::minstd_rand,
		class distType = std::uniform_real_distribution<inputType>,
		class resultType = inputType>
class RNG {
private:
	rngType engine;
	distType dist;
public:
	/*	Constructors */
	RNG(inputType first, inputType second, std::seed_seq s) :
			engine(s), dist(first, second) {
	}

	RNG(inputType first, inputType second,
			unsigned s =
					unsigned(
							std::chrono::system_clock::now().time_since_epoch().count())) :
			engine(s), dist(first, second) {
	}

	RNG(inputType param, std::seed_seq s) :
			engine(s), dist(param) {
	}

	RNG(inputType param,
			unsigned s =
					unsigned(
							std::chrono::system_clock::now().time_since_epoch().count())) :
			engine(s), dist(param) {
	}

	/*	These constructors are called when using discrete or peicewise distributions*/
	RNG(size_t n, double low, double high, double (*fw)(double),
			std::seed_seq s) :
			engine(s), dist(n, low, high, fw) {
	}

	RNG(size_t n, double low, double high, double (*fw)(double),
			unsigned s =
					unsigned(
							std::chrono::system_clock::now().time_since_epoch().count())) :
			engine(s), dist(n, low, high, fw) {
	}

	/*	Member functions	*/
	void set_seed(std::seed_seq s) {
		engine.seed(s);
	}
	void set_seed(unsigned s) {
		engine.seed(s);
	}

	void set_params(inputType first, inputType second = inputType()) {
		typename distType::param_type pt;
		if (second == inputType())
			pt = distType::param_type(first);
		else
			pt = distType::param_type(first, second);
		dist.param(pt);
	}

	/*	Called when using discrete or peicewise distributions */
	void set_params(size_t n, double xmin, double xmax, double (*fw)(double)) {
		typename distType::param_type pt(n, xmin, xmax, fw);
		dist.param(pt);
	}

	distType get_dist() const {
		return dist;
	}

	std::vector<resultType> random_vector(size_t n) {
		std::vector<resultType> retval;
		for (size_t i = 0; i < n; ++i)
			retval.push_back(dist(engine));
		return retval;
	}

	resultType random_number() {
		return dist(engine);
	}

};

// Due to ambiguities between constructor and member function signatures to other RNG types when templated,
// binomial distributions have been seperated out
template<class rngType = std::minstd_rand, class intType = int>
class RNG_binomial {
private:
	rngType engine;
	std::binomial_distribution<intType> dist;
public:
	RNG_binomial(intType upr_bnd, double prob, std::seed_seq s) :
			engine(s), dist(upr_bnd, prob) {
	}

	RNG_binomial(intType upr_bnd, double prob,
			unsigned s =
					unsigned(
							std::chrono::system_clock::now().time_since_epoch().count())) :
			engine(s), dist(upr_bnd, prob) {
	}

	void set_seed(std::seed_seq s) {
		engine.seed(s);
	}
	void set_seed(unsigned s) {
		engine.seed(s);
	}

	void set_params(intType upper, double prob) {
		typename std::binomial_distribution<intType>::param_type pt(upper,
				prob);
		dist.param(pt);
	}

	std::vector<intType> random_vector(size_t n) {
		std::vector<intType> retval;
		for (size_t i = 0; i < n; ++i)
			retval.push_back(dist(engine));
		return retval;
	}

	intType random_number() {
		return dist(engine);
	}
};

template<class rngType = std::minstd_rand, class intType = int>
class RNG_nbinomial {
private:
	rngType engine;
	std::negative_binomial_distribution<intType> dist;
public:
	RNG_nbinomial(intType upr_bnd, double prob, std::seed_seq s) :
			engine(s), dist(upr_bnd, prob) {
	}

	RNG_nbinomial(intType upr_bnd, double prob,
			unsigned s =
					unsigned(
							std::chrono::system_clock::now().time_since_epoch().count())) :
			engine(s), dist(upr_bnd, prob) {
	}

	void set_seed(std::seed_seq s) {
		engine.seed(s);
	}
	void set_seed(unsigned s) {
		engine.seed(s);
	}

	void set_params(intType upper, double prob) {
		typename std::negative_binomial_distribution<intType>::param_type pt(
				upper, prob);
		dist.param(pt);
	}

	std::vector<intType> random_vector(size_t n) {
		std::vector<intType> retval;
		for (size_t i = 0; i < n; ++i)
			retval.push_back(dist(engine));
		return retval;
	}

	intType random_number() {
		return dist(engine);
	}
};

} //namespace

#endif //header guard
