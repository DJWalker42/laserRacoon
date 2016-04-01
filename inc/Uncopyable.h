#ifndef UNCOPYABLE_HPP
#define UNCOPYABLE_HPP

class Uncopyable{
protected:
	Uncopyable(){}
	~Uncopyable(){}
private:
	Uncopyable( const Uncopyable& );
	Uncopyable& operator=( const Uncopyable& );
};

#endif
