#pragma once 

#include <array>

#define NUMCOUNTERS 5

namespace fem 
{

enum HPMType {
	vl_sum,
	v_instr, 
	cycles,
	fmem, 
	flops
}; 

/// class for reading hardware counters 
class HWCounter {
public:
	/// initialize for a given CSR 
	HWCounter() { }
	/// reset by storing current value 
	void Reset(); 
	/// get current values 
	void Read(); 

	/// return average vector length 
	double AvgVecLen() const; 
	/// return elapsed cycles 
	int Cycles() const; 
	/// return arithmetic intensity 
	double GetQ() const; 
	/// return flops per cycle 
	double FlopsPerCycle() const; 

	void PrintStats(std::string name="main") const; 
private:
	std::array<uint64_t,NUMCOUNTERS> _counters; 
}; 

} // end namespace fem 