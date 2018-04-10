#pragma once
/*
Particle - data structure to store position and velocity for one single particle
*/

template <class DataType>
class Particle
{
public:
	
	/* variable list */
	DataType x1;              // x1 - position
	DataType v1;              // v1 - velocity

	/* function list */
	Particle() :x1(0), v1(0) {}                          // constructor function 
	Particle(DataType x, DataType v) :x1(x), v1(v) {}    // alternative constructor function
};

/* sort funtion for Particle based on the first variable */
template <class DataType>
static bool particleSort(const Particle<DataType>& a1, const Particle<DataType>& a2)
{
	return a1.x1 <= a2.x1;
}