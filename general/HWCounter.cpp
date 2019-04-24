#include "HWCounter.hpp"
#include "General.hpp"

using namespace std; 

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
	if (_counters[HPMType::v_instr] == 0) return 0; 
	return (double)_counters[HPMType::vl_sum]/_counters[HPMType::v_instr];
}

int HWCounter::Cycles() const {
	return _counters[HPMType::cycles]; 
}

double HWCounter::GetQ() const {
	return (double)_counters[HPMType::flops]/_counters[HPMType::fmem]; 
}

double HWCounter::FlopsPerCycle() const {
	return (double)_counters[HPMType::flops]/_counters[HPMType::cycles]; 
}

void HWCounter::PrintStats(string name) const {
	cout << name << " stats:" << endl; 
	cout << "\taverage vl = " << AvgVecLen() << endl; 
	cout << "\tq = " << GetQ() << endl; 
	cout << "\tflops / cycle = " << FlopsPerCycle() << endl; 
}

} // end namespace fem 