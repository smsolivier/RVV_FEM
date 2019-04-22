.text
.align 2

#include "rvv.h"

.globl VectorScale_RV
.type  VectorScale_RV,@function

# N(0), alpha(1), v(2) 

VectorScale_RV:
	setvcfg(vcfg0,
      SCALAR | INT | W8,
      VECTOR | INT | W8,
		VECTOR | INT | W8,
      SCALAR | FP  | W64)
    setvcfg(vcfg2, \
      VECTOR | FP | W64, \
      VECTOR | FP | W64, \
      0, 0)
    vlds v3, a1 
loop: 
	setvl(t0, a0)
	vld v4, 0(a2)
	vmul v5, v4, v3 
	vst v5, 0(a2) 
	slli t1, t0, 3 
	add a2, a2, t1 
	sub a0, a0, t0 
	bnez a0, loop 
ret
