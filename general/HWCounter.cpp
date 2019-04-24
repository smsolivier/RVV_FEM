#include "HWCounter.hpp"
#include "General.hpp"

#define READ_ALL() \
	READ_CSR(hpmcounter3); \
	READ_CSR(hpmcounter4); \
	READ_CSR(cycle);\
	READ_CSR(hpmcounter5); \
	READ_CSR(hpmcounter6); 

namespace fem
{

void HWCounter::Reset() {
	int i=0; 
	#define READ_CSR(name) {\
		uint64_t val = read_csr(name); \
		_counters[i] = val; \
		i++; \
	}

	READ_ALL(); 
	#undef READ_CSR
}

void HWCounter::Read() {
	int i=0; 
	#define READ_CSR(name) {\
		uint64_t val = read_csr(name); \
		_counters[i] = val - _counters[i]; \
		i++; \
	}

	READ_ALL(); 
	#undef READ_CSR
}

double HWCounter::AvgVecLen() const {
	if (_counters[1] == 0) return 0; 
	return (double)_counters[0]/_counters[1];
}

int HWCounter::Cycles() const {
	return _counters[2]; 
}

double HWCounter::GetQ() const {
	return (double)_counters[3]/_counters[4]; 
}

} // end namespace fem 