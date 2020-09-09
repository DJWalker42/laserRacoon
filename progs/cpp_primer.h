#ifndef LASER_RACOON_EXAMPLE_CPP_PRIMER_H
#define LASER_RACOON_EXAMPLE_CPP_PRIMER_H

#include <string>
#include <iostream>

namespace example {

enum class season {winter, spring, summer, autumn};

/*toUType - function to return the underlying type of an enumerator
 * Typically used with scoped (class) enums but can be used with any enum structure
 */
template<typename E>
constexpr typename std::underlying_type<E>::type
toUType(E enumerator) noexcept{
  return static_cast<typename std::underlying_type<E>::type> (enumerator);
}

// the following structs and template functions fully implement the std::make_unique
// function found in C++14 upwards.
template<class T> struct _Unique_if {
    typedef std::unique_ptr<T> _Single_object;
};

template<class T> struct _Unique_if<T[]> {
    typedef std::unique_ptr<T[]> _Unknown_bound;
};

template<class T, size_t N> struct _Unique_if<T[N]> {
    typedef void _Known_bound;
};

template<class T, class... Args>
    typename _Unique_if<T>::_Single_object
    make_unique(Args&&... args) {
        return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
    }

template<class T>
    typename _Unique_if<T>::_Unknown_bound
    make_unique(size_t n) {
        typedef typename std::remove_extent<T>::type U;
        return std::unique_ptr<T>(new U[n]());
    }

template<class T, class... Args>
    typename _Unique_if<T>::_Known_bound
    make_unique(Args&&...) = delete;
// ----------------------------------------------------------------------------


//Abstract base class
class Item {
	//it is a preference to put internal representation (data members) first...
private:
	std::string _name;
	double _weight;
	int _quality;
	//.. constructors and destructor second ...
protected:
	//we cannot directly call this from code - must be done by derived classes
	Item(std::string name, double weight, int quality);

public:
	virtual ~Item(){} //need a virtual destructor as Item has a pure virtual function
	//... public interface third
public:
	//pure virtual function
	virtual void ability(season) = 0;

	//implicitly inlined member functions
	std::string name() const noexcept { return _name;}
	int quality() const noexcept {return _quality;}
	void setQuality(int new_quality) noexcept {_quality = new_quality;}

	//overload the output stream operator for an Item. Note this is a friend function
	//not a member function of Item. We return a ostream reference such that the operator
	//can be chained e.g. std::cout << an_item << another_item << std::endl;
	friend std::ostream& operator<<(std::ostream& os, const Item& item);

	//... and finally private internal functions, if any.
};

//Derived classes: Potion is an Item (but is still abstract, protected constructor)
class Potion : public Item {
protected:
	Potion(std::string name, int quality);
};

// Tonic is a Potion (concrete class, final prevents Tonic being a base class)
class Tonic final: public Potion {
private:
	//static data members are shared between objects of the same class
	static int _tonic_count; //initialised explicitly in .cpp file
public:
	//explicit to avoid implicit casts of primitive type int to a Tonic
	explicit Tonic(int quality);

public:
	void ability(season) override; //overriding the pure virtual function in Item

	//Note override can only be used for functions declared virtual, pure or otherwise.

	Tonic& operator+=(const Tonic& rhs);

	friend Tonic operator+(Tonic lhs, const Tonic& rhs);

	static int numberOfTonics() noexcept {return _tonic_count;}
};

// Poison is a Potion
class Poison final: public Potion {
public:
	explicit Poison(int quality);

public:
	void ability(season) override;
};


class Weapon : public Item {
protected:
	Weapon(std::string name, double weight, int quality);
};

//Sword is not "final" so could be used as a base class
class Sword : public Weapon {
public:
	Sword(double weight, int quality);

public:
	void ability(season) override;
};



}; //namespace

#endif //header guard
