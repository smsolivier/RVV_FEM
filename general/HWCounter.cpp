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
		_ctrs[i] = val; \
		i++; \
	}

	READ_ALL(); 
	#undef READ_CSR
}

void HWCounter::Read() {
	int i=0; 
	#define READ_CSR(name) {\
		uint64_t val = read_csr(name); \
		_ctrs[i] = val - _ctrs[i]; \
		i++; \
	}

	READ_ALL(); 
	#undef READ_CSR
}

double HWCounter::AvgVecLen() const {
	if (_ctrs[HPM::v_instr] == 0) return 0; 
	return (double)_ctrs[HPM::vl_sum]/_ctrs[HPM::v_instr];
}

int HWCounter::Cycles() const {
	return _ctrs[HPM::cycles]; 
}

double HWCounter::GetQ() const {
	if (_ctrs[HPM::fmem]==0) WARNING("memory access is zero"); 
	return (double)_ctrs[HPM::flops]/_ctrs[HPM::fmem]; 
}

double HWCounter::FlopsPerCycle() const {
	return (double)_ctrs[HPM::flops]/_ctrs[HPM::cycles]; 
}

uint64_t HWCounter::CacheAccesses() const {
	return _ctrs[HPM::accesses]; 
}

uint64_t HWCounter::CacheMisses() const {
	return _ctrs[HPM::misses]; 
}

double HWCounter::CacheMissRate() const {
	return (double)CacheMisses() / CacheAccesses() * 100.; 
}

uint64_t HWCounter::CacheBytesRead() const {
	return _ctrs[HPM::bytes_read]; 
}

double HWCounter::GetFOM(int tf, int lm, int lh) const {
	double m = CacheMissRate()/100; 
	return Cycles() + (lm*m + lh*(1-m))*CacheAccesses(); // + tf*_ctrs[HPM::flops]; 
}

void HWCounter::PrintStats(string name) const {
	cout << name << " stats:" << endl; 
	cout << "\taverage vl = " << AvgVecLen() << endl; 
	cout << "\tq = " << GetQ() << endl; 
	cout << "\tflops / cycle = " << FlopsPerCycle() << endl; 
	cout << "\tcycles = " << Cycles() << endl; 
	cout << "\tCache: " << CacheMisses() << "/" << CacheAccesses() << ", "
		<< (double)CacheMisses()/CacheAccesses() << "%, "
		<< CacheBytesRead() << "B read" << endl; 
	cout << "\tFOM = " << GetFOM() << " cycles" << endl; 
}

} // end namespace fem 