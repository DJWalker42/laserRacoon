/*
 *  A source file to highlight some of the features of the C++ language
 *  and STL
 *
 *  compile with: 	g++ -std=c++11 cpp_primer.cpp -o cpp_primer
 *  run with:		./cpp_primer [0|1|2|3]
 */
#include "cpp_primer.h"
#include <sstream>
#include <vector>
#include <map>
#include <stdexcept>
#include <algorithm>


namespace example {

//a look-up table for the season names, must match 'season' enumeration order and
// winter must be item zero
static const std::string season_name[] { "winter", "spring", "summer", "autumn"};

//a map might be better, we don't need to worry about order, the pairs are immutable, and
//we don't need to worry about the numbering of the enumeration
static const std::map<season, std::string> season_name_map {
		{season::winter, "winter"},
		{season::spring, "spring"},
		{season::autumn, "autumn"},
		{season::summer, "summer"}
};

Item::Item(std::string name, double weight, int quality) :
		_name(name), _weight(weight), _quality(quality)
{}

//overloaded operator not a member of Item so does not require the class name scope
std::ostream& operator<<(std::ostream& os, const Item& item) {
	os << item._name << " " << item._weight << " " << item._quality << "\n";
	return os;
}

Potion::Potion(std::string name, int quality) :
		Item(name, 0.1, quality)
{}

int Tonic::_tonic_count = 0;

Tonic::Tonic(int quality) : Potion("Tonic", quality) {
	_tonic_count++;
}

void Tonic::ability(season time_of_year) {
	int modifier {0};
	//use the LUT for the season names
	std::string the_season = season_name[toUType(time_of_year)];

	switch (time_of_year) {
	case season::winter:
		modifier = 1;
		break;
	case season::spring: //fall through wanted
	case season::summer: //summer and spring share the same modifier
		modifier = 3;
		break;
	case season::autumn:
		modifier = 0;
		break;
		//default unnecessary; a "season" can only be one of the specified four
	}

	std::string low_high = this->quality() < 75 ? "low" : "high";

	std::cout << "It is " << the_season << ": " << low_high << " quality "
			<< this->name() << " heals for " << modifier * this->quality()
			<< " hit points" << std::endl;
}

Tonic& Tonic::operator+=(const Tonic& rhs) {
	//data members in base class Item have private access,
	//so have to use public getter and setter functions
	this->setQuality(this->quality() + rhs.quality());
	return *this;
}

//passing lhs by value here helps optimise chained addition
Tonic operator+(Tonic lhs, const Tonic& rhs) {
	lhs += rhs; //reuse compound increment operator (above)
	return lhs; //return result by value (uses move constructor)
}

Poison::Poison(int quality) : Potion("Poison", quality) {}

void Poison::ability(season the_season) {
	if (the_season == season::summer) {
		throw std::runtime_error("Poison spoils in the heat");
	} else {
		std::cout << this->name() << " damages you for "
				<< this->quality() << " hit points" << std::endl;
	}
}

Weapon::Weapon(std::string name, double weight, int quality) :
		Item(name, weight, quality)
{}

Sword::Sword(double weight, int quality) : Weapon("Sword", weight, quality) {}

void Sword::ability(season time_of_year) {

	//use the map to look-up the season name
	std::string the_season = season_name_map.find(time_of_year)->second;

	std::ostringstream os;

	switch (time_of_year) {
	case season::winter: //fall through wanted
	case season::spring: //winter and spring share the same effect
		os << "It is " << the_season << " your " << this->name()
		<< " is frozen in its scabbard - run away!\n";
		break;
	case season::summer:
		os << "It is " << the_season
		<< " your hand is sweaty and the hilt slips from your grip - run away!\n";
		break;
	case season::autumn:
		os << "It is " << the_season << " your " <<
		this->name() << " has rusted in its scabbard from lack of use - run away!\n";
		break;
		//default unnecessary; a "season" can only be one of the specified four
	}

	std::cout << os.str() << std::endl;
}


} //namespace


//char ** argv is equivalent to char * argv[]
int main (int argc, char ** argv) {

	using namespace example; //because I'm a lazy typer

	using UPtrItem = std::unique_ptr<Item>;
	using Backpack = std::vector<UPtrItem>;

	int choice_of_season = 0;

	const int q_bound = 75; //item quality boundary

	if (argc > 1) {
		choice_of_season = atoi(argv[1]);
	}

	if (choice_of_season > 3) {
		std::cerr << "There are only 4 seasons! Using winter\n";
		choice_of_season = 0;
	}

	if (choice_of_season < 0) {
		std::cerr << "Negative season does not compute! Using summer\n";
		choice_of_season = 2;
	}

	//there is no direct conversion from an integer to a enumeration type
	season the_season;
	switch(choice_of_season) {
	case 0:
		the_season = season::winter;
		break;
	case 1:
		the_season = season::spring;
		break;
	case 2:
		the_season = season::summer;
		break;
	case 3:
		the_season = season::autumn;
		break;
		//no default - already dealt with invalid choices
	}

	//we can call a static member function without the need to instantiate a Tonic object
	std::cout << "Number of Tonics created: " << Tonic::numberOfTonics() << std::endl;

	Tonic low_quality(10);
	Tonic high_quality(100);

	//because we have overload the addition operator for Tonic we can write the following:
	Tonic mix = low_quality + high_quality;

	//Unless we add similar overloads to other derived classes, only Tonic objects can be added

	//also, as we overloaded the output stream operator for Item base class we can write:
	std::cout << "Mixed Tonic: " << mix << std::endl;

	//Will it be 2 or 3 here?
	//Note that the addition assignment uses the compiler provided copy constructor
	std::cout << "Number of Tonics created: " << Tonic::numberOfTonics() << std::endl;

	//hero gets an empty backpack
	Backpack backpack;

	//hero picks up some items
	backpack.push_back(make_unique<Tonic>(100));
	backpack.push_back(make_unique<Tonic>(50));
	backpack.push_back(make_unique<Poison>(20));
	backpack.push_back(make_unique<Sword>(8.7, 1000));

	//output each item placed in the backpack
	for(auto& item : backpack) {
		std::cout << *item; //dereference the unique_ptr to get an Item object
	}
	std::cout << std::endl;

	//has the hero picked up any "low quality" items?
	//Here we use a lambda to express the predicate for the algorithm 'std::any_of'
	if (std::any_of(backpack.begin(), backpack.end(),
			[q_bound](const UPtrItem& item){return item->quality() < q_bound;})) {
		std::cout << "Backpack contains at least one low quality item" << std::endl;
	}

	//hero uses each of the items in their backpack
	for(auto& item : backpack) {
		//try-catch block for exceptions
		try {
			//item is a unique_ptr to an Item object
			//there is no consideration of the specific derived object
			//we just call the interface function for each item
			item->ability(the_season);
		} catch (const std::runtime_error& e){
			std::cerr << "Something exceptional: " << e.what() << std::endl;
			//code continues from after the catch block when an exception is thrown
		} catch (...) {
			//the ellipsis ... means catch anything else not covered by runtime_error
		}
	}

	std::cout << "Number of Tonics created: " << Tonic::numberOfTonics() << std::endl;

	return 0;
}






