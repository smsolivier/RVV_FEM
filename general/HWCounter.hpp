#pragma once 

#include <array>

#define NUMCOUNTERS 8

namespace fem 
{

/// keeps track of ordering of the _ctrs array 
enum HPM {
	vl_sum,
	v_instr, 
	cycles,
	fmem, 
	flops,
	accesses, 
	misses, 
	bytes_read
}; 

/// class for reading hardware counters 
class HWCounter {
public:
	/// initialize for a given CSR 
	HWCounter() {Reset(); }
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
	/// return cache accesses 
	uint64_t CacheAccesses() const; 
	/// return cache misses 
	uint64_t CacheMisses() const; 
	/// return bytes read 
	uint64_t CacheBytesRead() const; 
	/// return the cache miss rate (percentage) 
	double CacheMissRate() const; 

	/// return an estimate of execution time 
	/** \param tf number of cycles per flop 
		\param lm latency of cache miss 
		\param lh latency of cache hit 
	*/ 
	double GetFOM(int tf=3, int lm=20, int lh=1) const; 

	void PrintStats(std::string name="main") const; 
private:
	std::array<uint64_t,NUMCOUNTERS> _ctrs; 
}; 

} // end namespace fem 