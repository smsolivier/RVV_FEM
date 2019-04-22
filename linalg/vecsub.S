.text
.align 2

#include "rvv.h"

.globl VectorSub_RV
.type  VectorSub_RV,@function

# N(0), a(1), b(2) 

VectorSub_RV:
	setvcfg(vcfg0, 
		VECTOR | FP | W64, 
		VECTOR | FP | W64,
		VECTOR | FP | W64,
		VECTOR | FP | W64
	)
loop:
	setvl(t0, a0) 
	vld v0, 0(a1) 
	vld v1, 0(a2) 
	vsub v0, v0, v1 
	vst v0, 0(a1) 
	slli t1, t0, 3 
	add a1, a1, t1 
	add a2, a2, t1 
	sub a0, a0, t0 
	bnez a0, loop 
ret 
