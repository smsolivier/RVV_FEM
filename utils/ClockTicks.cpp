#include "ClockTicks.hpp"
#include "General.hpp"

#ifdef __riscv 
uint64_t CH_pticks; 
unsigned long long int ch_ticks()
{
  uint64_t ticks = read_csr(cycle); 
  uint64_t diff = ticks - CH_pticks; 
  CH_pticks = ticks; 
  // return diff;
  return ticks;  
}

uint64_t CH_vl = 0; 
uint64_t CH_vi = 0; 
uint64_t CH_flops = 0; 
uint64_t CH_mem = 0; 
double ch_avl() {
	uint64_t vl = read_csr(hpmcounter3); 
	uint64_t vi = read_csr(hpmcounter4); 
	uint64_t vle = vl - CH_vl; 
	uint64_t vie = vi - CH_vi; 
	CH_vl = vl; 
	CH_vi = vi; 

	if (vie == 0) return 0; 
	return (double)vle/vie; 
}

double ch_q() {
	uint64_t mem = read_csr(hpmcounter5); 
	uint64_t flops = read_csr(hpmcounter6); 
	uint64_t meme = mem - CH_mem; 
	uint64_t flopse = flops - CH_flops; 
	CH_mem = mem; 
	CH_flops = flops; 
	if (flopse==0 or meme==0) return 0; 
	return (double)flopse/meme; 
}
#else 
double ch_avl() {return 0; }
double ch_q() {return 0.; }
#endif
