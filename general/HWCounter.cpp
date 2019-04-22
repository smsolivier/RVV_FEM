#include "HWCounter.hpp"

namespace fem
{

void HWCounter::Reset() {
	int i=0; 
	#define READ_CSR(name) {\
		uint64_t val = read_csr(name); \
		_counters[i] = val; \
		i++; \
	}

	READ_CSR(hpmcounter3); 
	READ_CSR(hpmcounter4); 
	#undef READ_CSR
}

void HWCounter::Read() {
	int i=0; 
	#define READ_CSR(name) {\
		uint64_t val = read_csr(name); \
		_counters[i] = val - _counters[i]; \
		i++; \
	}

	READ_CSR(hpmcounter3); 
	READ_CSR(hpmcounter4); 
	#undef READ_CSR
}

double HWCounter::AvgVecLen() const {
	return (double)_counters[0]/_counters[1];
}

} // end namespace fem 