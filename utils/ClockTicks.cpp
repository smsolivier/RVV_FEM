#include "ClockTicks.hpp"
#include "General.hpp"

uint64_t CH_pticks; 
unsigned long long int ch_ticks()
{
  uint64_t ticks = read_csr(cycle); 
  uint64_t diff = ticks - CH_pticks; 
  CH_pticks = ticks; 
  // return diff;
  return ticks;  
}

uint64_t CH_vl; 
uint64_t CH_vi; 
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
