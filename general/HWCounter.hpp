#pragma once 

#include "Array.hpp"
#include <array>

#define NUMCOUNTERS 3

namespace fem 
{

/// class for reading hardware counters 
class HWCounter {
public:
	/// initialize for a given CSR 
	HWCounter() { }
	/// reset by storing current value 
	void Reset(); 
	/// get current values 
	void Read(); 

	/// get average vector length 
	double AvgVecLen() const; 
	/// get elapsed cycles 
	int Cycles() const; 
private:
	std::array<uint64_t,NUMCOUNTERS> _counters; 
}; 

} // end namespace fem 