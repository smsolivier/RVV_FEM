#include "HWCounter.hpp"
#include "General.hpp"

using namespace std; 

#ifdef __riscv
#define read_csr(reg) ({ unsigned long __tmp; \
  asm volatile ("csrr %0, " #reg : "=r"(__tmp)); \
  __tmp; })
#else 
#define read_csr(reg) 0; 
#endif 

#define READ_ALL() \
	READ_CSR(hpmcounter3); \
	READ_CSR(hpmcounter4); \
	READ_CSR(cycle);\
	READ_CSR(hpmcounter5); \
	READ_CSR(hpmcounter6); \
	READ_CSR(hpmcounter7); \
	READ_CSR(hpmcounter8); \
	READ_CSR(hpmcounter9); 

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
	if (_counters[HPM::v_instr] == 0) return 0; 
	return (double)_counters[HPM::vl_sum]/_counters[HPM::v_instr];
}

int HWCounter::Cycles() const {
	return _counters[HPM::cycles]; 
}

double HWCounter::GetQ() const {
	if (_counters[HPM::fmem]==0) WARNING("memory access is zero"); 
	return (double)_counters[HPM::flops]/_counters[HPM::fmem]; 
}

double HWCounter::FlopsPerCycle() const {
	return (double)_counters[HPM::flops]/_counters[HPM::cycles]; 
}

uint64_t HWCounter::CacheAccesses() const {
	return _counters[HPM::accesses]; 
}

uint64_t HWCounter::CacheMisses() const {
	return _counters[HPM::misses]; 
}

uint64_t HWCounter::CacheBytesRead() const {
	return _counters[HPM::bytes_read]; 
}

void HWCounter::PrintStats(string name) const {
	cout << name << " stats:" << endl; 
	cout << "\taverage vl = " << AvgVecLen() << endl; 
	cout << "\tq = " << GetQ() << endl; 
	cout << "\tflops / cycle = " << FlopsPerCycle() << endl; 
	cout << "\tCache: " << CacheMisses() << "/" << CacheAccesses() << ", "
		<< (double)CacheMisses()/CacheAccesses() << "%, "
		<< CacheBytesRead() << "B read" << endl; 
}

} // end namespace fem 